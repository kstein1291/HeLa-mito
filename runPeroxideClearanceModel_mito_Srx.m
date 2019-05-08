function [out,C] = runPeroxideClearanceModel_mito_Srx
%baseline model for hydrogen peroxide reaction network in the mitochondria
%Kassi T. Stein
%June 2018

%The overall function performs a parameter sweep over possible
%concentrations of Prx3 and solves for the baseline concentrations of 28
%different mitochondrial species, based on a constant H2O2 source from
%OxPhos. The function outputs are "C" a vector of Prx3 concentrations from
%the parameter sweep and "out" a cell array containing a matrix of the
%output from each iteration of the function in each cell.

out=cell(1,10);
C=linspace(48,110,10);
for m=1:length(out)
    [OUT]=PeroxideClearanceModel_HeLa_mito(C(m),m);
    out{m}=OUT;
end

end

function out = PeroxideClearanceModel_HeLa_mito(C,q)
% Model for baseline hydrogen peroxide consumption
% Kassi T. Stein
% June 2018
% model takes as input a concentration of Prx3 and an iteration number (for
% making the Excel file output)
% model outputs a matrix "out" where the first column is time and each 
% susequent column is each component of the vector x (x1, x2,...,x28)
% the model also outputs and Excel file where each iteration of the
% parameter sweep is its own worksheet (make sure to specify a filename)

% Define model parameters
k = zeros(29,1);
% Intracellular peroxide production
k(1) = 4; %uM/s OCR est
% GPx1red reacting with H2O2
k(2) = 60; % uM^-1*s^-1
% GPx1ox reacting with GSH
k(3) = 4e-2; % uM^-1*s^-1
% GPx-SSG reacting with GSH
k(4) = 10; % uM^-1*s^-1 (assuming same for both isoforms)
% Km of NADP+
k(5) = 57; % uM
% Prx3-SH oxidized by H2O2
k(6) = 20; % uM^-1*s^-1
% Prx3-SOH over-oxidized by H2O2
k(7) = 1.4e-2; % uM^-1*s^-1
% Reduction of overoxidized Prx3 by Srx enzyme
k(8) = 3e-3; % uM^-1*s^-1 
% Self-catalyzed disulfide formation of Prx3-SS from Prx3-SOH
k(9) = 20; % s^-1
% Prx3 is reduced by thioredoxin
k(10) = 2.2e-1; % uM^-1*s^-1
% Auto-oxidation of GSH
k(11) = 7.4e-05; % s^-1
% Pr-SH oxidized by H2O2
k(12) = 1e-4; % uM^-1*s^-1
% Pr-SOH glutathionylated by GSH
k(13) = 1.2e-1; % uM^-1*s^-1
% Grx2-SH de-glutathionylates Protein-SSG
k(14) = 1.2e-2; % uM^-1*s^-1
% GSH de-glutathionylates Grx2-SSG
k(15) = 3.7e-2; % uM^-1*s^-1
% Pr-(SH)2 oxidized by H2O2
k(16) = 1e-4; % uM^-1*s^-1
% Pr-SS reduced by Trx
k(17) = 1e-4; % uM^-1*s^-1
% GSSG reduced by GR
k(18) = 3.2; % uM^-1*s^-1
% Oxidized Thioredoxin reduced by TrxR
k(19) = 20; % uM^-1*s^-1
% Production of NADPH by G6P-DH
k(20) = 3.75e2; % uM/s 
% GSH synthesis
k(21) = 4.1e-1; % uM/s
% Trx_SH synthesis
k(22) = 6.97e-4; % uM/s
%Prx5-SH oxidized by H2O2
k(23) = 3e-1; %uM^-1 s^-1
%Prx5-SOH auto-catalyzes to Prx5-SS
k(24) = 14.7; %s^-1
%Prx5-SS reduced by Trx2
k(25) = 2; % uM^-1 s^-1
%GPx4red oxidized by H2O2
k(26) = 4.8e-2; % uM^-1 s^-1
%Gpx4ox reacting with GSH
k(27) = 2e-2; % uM^-1 s^-1
%DAAO producing H2O2
k(28) = 0; %uM/s (we'll use this in the perturbation case)
%Srx import
k(29) = 1.23e-5; %uM/s

%%

% Defining initial Conditions (Concentrations in uM = micromoles/Liter)
x0 = zeros(28,1);
x0(1) = .484; %H2O2 (this is a placeholder, we'll solve for it later)
x0(2) = 1.5e-2; %Gpx1red
x0(5) = 5e3; % GSH
x0(6) = 1.78; % GSSG 
x0(7) = C; % Prx3-SH
x0(11) = 7.7; % Trx2-SH
x0(12) = 7.54e-2; % Trx2SS
x0(13) = 1e-3; % Pr-SH 
x0(14) = x0(13)*(.5/100); % Pr-SOH
x0(15) = x0(13)*(.5/100); % Pr-SSG
x0(16) = 1; %Grx2-SH
x0(17) = x0(16)*(.5/100); % Grx-SSG
x0(18) = 1.09e3; % Pr-(SH)2
x0(20) = 30; % NADPH
x0(21) = 3.0e-1; % NADP+
x0(22) = 14; %Prx5-SH
x0(25) = 0.23; %Gpx4red
x0(28) = 8.78e-3; %Srx

% Specify baseline H2O2 concentration

initial_H2O2 = k(1)/(k(6)*x0(7));
x0(1) = initial_H2O2;
% Solve for expected initial conditions for GPXo and GPX-SG
x0(3) = x0(2)/(k(3)*x0(5)/k(2)/initial_H2O2+1+k(2)/k(4)); % GPXo
x0(4) = x0(2)/(k(4)*x0(5)/k(2)/initial_H2O2+k(4)/k(3)+1); % GPX-SG
x0(2) = x0(2) - x0(3) - x0(4); % Resolving for GPXr based on molar balance

% Solve for expected initial conditions for different oxidation states
% of Prx3
x0(8) = x0(7)/(k(9)/k(6)/initial_H2O2+1+k(9)/k(10)/x0(11)+k(7)*initial_H2O2/k(8)/x0(28));% Prx-SOH
x0(9) = x0(7)*k(7)*initial_H2O2/(k(8)*x0(28)*(k(9)/k(6)/initial_H2O2+1+k(9)/k(10)/x0(11)+k(7)*initial_H2O2/k(8)/x0(28))); % Prx-SOOH
x0(10) = x0(7)*k(9)/(k(10)*x0(11)*(k(9)/k(6)/initial_H2O2+1+k(9)/k(10)/x0(11)+k(7)*initial_H2O2/k(8)/x0(28))); % Prx-SS
x0(7) = x0(7) - x0(8) - x0(9) - x0(10); % Resolving for Prx-SH based on molar balance


% Solve for expected initial conditions for different oxidation states
% of Pr-SH
x0(14) = x0(13)/(k(13)*x0(5)/k(12)/initial_H2O2+k(13)*x0(5)/k(14)/x0(16)+1); % Pr-SOH
x0(15) = x0(13)/(k(14)*x0(16)/k(12)/initial_H2O2+1+k(14)*x0(16)/k(13)/x0(5)); % Pr-SSG
x0(13) = x0(13) - x0(14) - x0(15); % Resolving for Pr-SH based on molar balance

% Solve for expected initial condition for oxidized Pr-(SH)2, then use
% molar balance to resolve for Pr-(SH)2
x0(19) = x0(18)/(1+k(17)*x0(11)/k(16)/initial_H2O2); % Pr-SS
x0(18) = x0(18) - x0(19);

% Solve for expected initial condition for oxidized Grx-SSG, then use
% molar balance to resolve for Grx-SH
x0(17) = k(16)/(1+k(15)*x0(5)/k(14)/x0(15));
x0(16) = x0(16) - x0(17);

%solve for expected initial conditions of Prx5 states, then use molar
%balance to find Prx5-SH remaining
x0(23) = k(23)*x0(22)*initial_H2O2/k(24); %Prx5-SOH
x0(24) = k(23)*x0(22)*initial_H2O2/(k(25)*x0(11)); %Prx5-SS
x0(22) = x0(22)-x0(23)-x0(24); %Prx5-SH

% Solve for expected initial conditions for GPX4o and GPX4-SG
x0(26) = x0(26)/(k(27)*x0(5)/k(26)/initial_H2O2+1+k(26)/k(4)); % GPX4o
x0(27) = x0(26)/(k(4)*x0(5)/k(26)/initial_H2O2+k(4)/k(27)+1); % GPX4-SG
x0(25) = x0(25) - x0(26) - x0(27); % Resolving for GPX4r based on molar balance

%%
% Solver Parameters
ti = 0; % start time
tf = 5; % stop time (s)
%ode15s is a stiff ODE solver, so will use variable time steps at short vs.
%long times; only need to specify start and end time, not time steps
tspan = [ti tf];
%force the function to be always positive
negspan = [1:28];
options=odeset('AbsTol',1E-10,'RelTol',1E-4,'NonNegative',negspan);

% Integration
tic
[t,x]=ode15s(@crank,tspan,x0,options,k);
toc
out = [t, x];
filename = 'HeLa_srxn1.xlsx';
xlswrite(filename,out,q,'B3');

%%
% Description of derivatives
function dxdt = crank(t, x, k);
dxdt= x; % setting up vector containing derivatives
dxdt(1) = k(1) - k(2)*x(2)*x(1) - k(6)*x(7)*x(1) - k(7)*x(8)*x(1)... 
    - k(12)*x(13)*x(1) - k(16)*x(18)*x(1)-k(23)*x(22)*x(1)-k(26)*x(25)*x(1); % H2O2
dxdt(2) = -k(2)*x(2)*x(1) + k(4)*x(4)*x(5); % GPX1red
dxdt(3) = k(2)*x(2)*x(1) - k(3)*x(3)*x(5); % GPX1ox
dxdt(4) = k(3)*x(3)*x(5) - k(4)*x(4)*x(5); % GPX1-SG
dxdt(5) = -k(3)*x(3)*x(5) - k(4)*x(4)*x(5) - 2*k(11)*x(5)... 
    - k(13)*x(14)*x(5) - k(15)*x(17)*x(5) + 2*k(18)*x(6)*x(20) + k(21)...
    - k(27)*x(26)*x(5) - k(4)*x(27)*x(5); % GSH
dxdt(6) = k(4)*x(4)*x(5) + k(11)*x(5) + k(15)*x(17)*x(5) + k(4)*x(27)*x(5) - k(18)*x(6)*x(20); % GSSG
dxdt(7) = -k(6)*x(7)*x(1) + k(10)*x(10)*x(11); % Prx3-SH
dxdt(8) = k(6)*x(7)*x(1) - k(7)*x(8)*x(1) + k(8)*x(9)*x(28) - k(9)*x(8); % Prx-SOH
dxdt(9) = k(7)*x(8)*x(1) - k(8)*x(9)*x(28); % Prx-SOOH
dxdt(10) = k(9)*x(8) - k(10)*x(10)*x(11); % Prx3-SS
dxdt(11) = -k(10)*x(10)*x(11) - k(17)*x(19)*x(11) - k(25)*x(24)*x(11)...
    + k(19)*x(12)*x(20) + k(22); % Trx-SH
dxdt(12) = k(10)*x(10)*x(11) + k(17)*x(19)*x(11) + k(25)*x(24)*x(11)...
    - k(19)*x(12)*x(20); % Trx-SS
dxdt(13) = -k(12)*x(13)*x(1) + k(14)*x(16)*x(15); % Pr-SH
dxdt(14) = k(12)*x(13)*x(1) - k(13)*x(14)*x(5); % Pr-SOH
dxdt(15) = k(13)*x(14)*x(5) - k(14)*x(16)*x(15); % Pr-SSG
dxdt(16) = k(15)*x(17)*x(5) - k(14)*x(16)*x(15); % Grx-SH
dxdt(17) = k(14)*x(16)*x(15) - k(15)*x(17)*x(5); % Grx-SSG
dxdt(18) = -k(16)*x(18)*x(1) + k(17)*x(19)*x(11); % Pr-(SH)2
dxdt(19) = k(16)*x(18)*x(1) - k(17)*x(19)*x(11); % Pr-SS
dxdt(20) = -k(18)*x(6)*x(20) - k(19)*x(12)*x(20) + k(20)*x(21)/(k(5) + x(21)); % NADPH
dxdt(21) = k(18)*x(6)*x(20) + k(19)*x(12)*x(20) - k(20)*x(21)/(k(5) + x(21)); % NADP+
dxdt(22) = -k(23)*x(22)*x(1) + k(25)*x(24)*x(11); %Prx5-SH
dxdt(23) = k(23)*x(22)*x(1) - k(24)*x(23); %Prx5-SOH
dxdt(24) = k(24)*x(23) - k(25)*x(24)*x(11); %Prx5-SS
dxdt(25) = -k(26)*x(25)*x(1) + k(4)*x(27)*x(5); % GPX4red
dxdt(26) = k(26)*x(25)*x(1) - k(27)*x(26)*x(5); % GPX4ox
dxdt(27) = k(27)*x(26)*x(5) - k(4)*x(27)*x(5); % GPX4-SG
dxdt(28) = k(29); %Srx
end

end