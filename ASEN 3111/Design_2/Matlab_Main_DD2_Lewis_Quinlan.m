%% ASEN 3111 - Design Deliverable 2 - Main
%The purpose of this code is to predict a drag polar acting on a Tempest
%glider by using Theories used in class. This code uses all of the previous
%constants given from the last design deliverable as the ones given for
%this design deliverable. 

% Author: Quinlan Lewis
% Date: 11/5/20
% Last Modified: 11/5/20

%% Housekeeping
clear all; clc; close all

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

%% Computing CD_0 and k constant
%calculating k as a constant, N represents how many values wanted for the
%creating of the c array
k = 1/(pi*e*AR)
N = 100;

%from flat plat approximation we know that Re_x_cr = 5*10^5
Re_x_cr = 5 * 10^5;

%calculating all of the constants known as well as retrieving any values
%needed in these calculations from the function atmos that was found on the
%MATHWORKS website
V_inf = mean(linspace(17,25,N));
[rho_inf,a,T,P,mu_inf,z,sigma] = atmos(h);
q_inf = .5*rho_inf*V_inf^2;
x_cr = (mu_inf)*(Re_x_cr)/(V_inf*rho_inf);
[c] = ChordArrayFunc(1.53,N);
c = c(c~=0);
% x_cr = .25*(max(c))

%defines values for all of the Cf coefficients
Cf_1_lam = 1.328/sqrt(Re_x_cr);
Cf_1_turb = 0.074/((Re_x_cr)^(1/5));
Cf_c_turb = 0.074./((rho_inf*V_inf.*c./mu_inf).^(1/5));

%this below if statement is for the case when the transition point is
%greater than the chord length of the wings
if x_cr > max(c)
    %transition to turbulent flow happens off of wing so all of flow is
    %laminar
    Df_1_lam = q_inf*c*Cf_1_lam;
    Df_2_turb = 0;
    Df = Df_2_turb + Df_1_lam;
else
    %transition to turbulent flow happens on the wing so we need to take
    %into account the chord length and the transition point
    Df_2_turb = (q_inf.*c.*Cf_c_turb) - q_inf*x_cr*Cf_1_turb;
    Df_1_lam = q_inf*x_cr*Cf_1_lam;
    Df = Df_2_turb + Df_1_lam;
end

CD_0 = 2*mean(Df./(q_inf.*c))

%below are the calculations of gross weight to be used in the max range and
%max endurance calculations, as well as the velocity needed for max range
%and endurance
Wgross = W_0*9.81;
U_E = sqrt((2*Wgross)/(rho_inf*S)*(sqrt(k/(3*CD_0))));
U_R = sqrt((2*Wgross)/(rho_inf*S)*(sqrt(k/(CD_0))));

%Below are the equations used in calculating maximum range and endurance
E_max = ((Rt)^(1-n)) * ((n_tot*Volt*Cap) / ((2/sqrt(rho_inf*S)) * (CD_0^(1/4)) * (2*Wgross*sqrt(k/3))^(3/2)))^n
R_max = (((Rt)^(1-n)) * ((n_tot*Volt*Cap) / ((1/sqrt(rho_inf*S)) * (CD_0^(1/4)) * (2*Wgross*sqrt(k))^(3/2)))^n)*U_R*3.6