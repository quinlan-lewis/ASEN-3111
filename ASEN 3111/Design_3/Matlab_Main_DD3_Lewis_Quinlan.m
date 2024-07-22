%% ASEN 3111 - Design Deliverable 3 - Main
%The purpose of this code is to predict a drag polar acting on a Tempest
%glider by using Theories used in class. This code uses all of the previous
%constants given from the last design deliverable as the ones given for
%this design deliverable. 

% Author: Quinlan Lewis
% Date: 12/6/20
% Last Modified: 12/5/20

%% Housekeeping
clear; clc; close all;

%% Constants and Data
%from PLLT we can find that the oswalds effiecency number is e = 1 for an
%elliptical shaped wing
e = 1;

%setting given constants for geometry and physical characteristics
h = 1829; %cruising altitude in meters
n_tot = 0.5; %overall propulsion effeciency
b = 3.22; %wingspan
S = 0.63; %wing area
AR = 16.5; %aspect ratio
%airfoil = MH 32;
l = 1.56; %length
W = 4.5; %empty weight
W_0 = 6.4; %gross takeoff weight

%Constants for batter characteristics
%type = 2 packs, lithium polymer(Li-Po)
n = 1.3; %common discharge ratte for Li-Po battery
Cap = 18; %capacity [Ah]
Volt = 26.1; %voltage [V]
Rt = 1; %rating in h
BP_m = 2*(984/1000); %battery pack mass, 2 of them in grams
cur_max = 2*(50); %2 motors with 50 amps per motor

%% Computing CD_0 constant
%calculating k as a constant, N represents how many values wanted for the
%creating of the c array
%elliptical estimate for taper ratio = .35
N = 100;
V_inf = 27;
[rho_inf,a,T,P,mu_inf,z,sigma] = atmos(h);
MH32;
Data = airfoil.Body;
XB = Data(:,1);
YB = Data(:,2);
VINF = mean(V_inf);

for ALPHA = -5:1:5
    [CL(1,6+ALPHA),~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
end
a = (2*pi)/(1+2/AR);

ALPHA = -5:1:5;
Fit = polyfit(ALPHA'*pi/180,squeeze(CL(1,:)),1);

a0(1) = Fit(1);
alpha_0(1) = (-Fit(2)/Fit(1))*(180/pi);

b = b; %[ft]
c_r = .23; %[ft] root chord
c_t = c_r*.35; %[ft] tip chord
a0_r = a0(1); %[per rad]
a0_t = a0(1); %[per rad]
aero_r = (alpha_0(1)); %[rad] zero lift aoa
aero_t = (alpha_0(1)); %[rad] zero lift aoa
geo_r = deg2rad(mean(CL/a + alpha_0(1))); %[rad] geo aoa
geo_t = geo_r; %[rad] geo aoa
N = 100; % number of odd terms

%PLLT function is ran to get e and the corresponding coeffs
[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,0);
k(1) = 1/(pi*e*AR);
CD_0(1) = (c_L^2)/(pi*e*AR);

fprintf('geometric angle of attack for task 1: %f\n', geo_r);

%% Computing k constant using PLLT and Vortex Panel Method

%Obtaining a0 and alpha_L=0 for MH32
%CA3 Solution code is used to get values for VP_Slope(1) and
%VP_ZeroLiftAOA(1) for the tip values

MH32;
Data = airfoil.Body;
XB = Data(:,1);
YB = Data(:,2);
VINF = mean(V_inf);

for ALPHA = -5:1:5
    [CL(1,6+ALPHA),~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
end

ALPHA = -5:1:5;
Fit = polyfit(ALPHA'*pi/180,squeeze(CL(1,:)),1);

a0(1) = Fit(1);
alpha_0(1) = (-Fit(2)/Fit(1))*(180/pi);

%using above values I can use my PLLT function to calculate the c_L of the
%tempest

b = b; %[ft]
c_r = .23; %[ft] root chord
c_t = .23; %[ft] tip chord
a0_r = a0(1); %[per rad]
a0_t = a0(1); %[per rad]
aero_r = (alpha_0(1)); %[rad] zero lift aoa
aero_t = (alpha_0(1)); %[rad] zero lift aoa
geo_r = deg2rad(mean(CL/a0(1) + alpha_0(1))); %[rad] geo aoa
geo_t = geo_r; %[rad] geo aoa
N = 100; % number of odd terms

%PLLT function is ran to get e and the corresponding coeffs
[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,0);
k(2) = 1/(pi*e*AR);
CD_0(2) = (c_L^2)/(pi*e*AR);

fprintf('geometric angle of attack for task 2: %f\n', geo_r);
%% Calculating range and endurance
Wgross = W_0*9.81;
U_E = sqrt((2.*Wgross)./(rho_inf.*S).*(sqrt(k./(3.*CD_0))));
U_R = sqrt((2.*Wgross)./(rho_inf.*S).*(sqrt(k./(CD_0))));

%Below are the equations used in calculating maximum range and endurance
E_max = ((Rt).^(1-n)) .* ((n_tot.*Volt.*Cap) ./ ((2./sqrt(rho_inf.*S)) .* (CD_0.^(1/4)) .* (2.*Wgross.*sqrt(k./3)).^(3/2))).^n;
R_max = (((Rt).^(1-n)) .* ((n_tot.*Volt.*Cap) ./ ((1./sqrt(rho_inf.*S)) .* (CD_0.^(1/4)) .* (2.*Wgross.*sqrt(k)).^(3/2))).^n).*U_R.*3.6;

fprintf('Task 1 values:\n k = %f\n cD_0 = %f\n U_E = %f [m/s]\n U_R = %f [m/s]\n E_max = %f [hrs]\n R_max = %f [kilometers]\n',k(1), CD_0(1), U_E(1), U_R(1), E_max(1), R_max(1))
fprintf('Task 2 values:\n k = %f\n cD_0 = %f\n U_E = %f [m/s]\n U_R = %f [m/s]\n E_max = %f [hrs]\n R_max = %f [kilometers]\n',k(2), CD_0(2), U_E(2), U_R(2), E_max(2), R_max(2))