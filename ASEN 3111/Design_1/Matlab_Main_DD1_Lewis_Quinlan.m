%% ASEN 3111 - Design Deliverable 1 - Main
%Purpose of this code is to analyze given CFD data and and use that data to
%determine values for CD_0 and k, then using those determined values,
%determine the maximum range and endurance off of the given data and
%constants.

% Author: Quinlan Lewis
% Date: 10/15/20
% Last Modified: 10/15/20

%% Housekeeping
clear all; clc; close all

%% Constants and Data
CFD_data = readmatrix('3111 DLD#1 CFD - Sheet1.csv');
alpha = CFD_data(:,1);
CL = CFD_data(:,2);
CD = CFD_data(:,3);

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
%CD_0 occurs at zero lift which is where CL = 0, could use polyval and
%polyfit to determine a funciton then evaluate that function at CL = 0
%can use funciton atmoscoesa() to find temperature, pressure, speed of
%sound, and rho at a given altitube in meters

%Atmoscoesa gives the density, pressure, speed of sound, and temperature at
%a given altitude h
[T, a, P, rho] = atmoscoesa(h);

CL_trimmed = CL(4:14);
CD_trimmed = CD(4:14);
alpha_trimmed = alpha(4:14);
%Below I create different polyfits to create lines of best fit to be
%evaluated for CD_0 and k
p_CLvCD = polyfit(CL_trimmed,CD_trimmed,2);
p_CL_2vCD = polyfit(CL_trimmed.^2,CD_trimmed,1);
CLvAlpha = polyfit(alpha_trimmed,CL_trimmed,1);
CDvAlpha = polyfit(alpha_trimmed,CD_trimmed,2);

%Below are 4 figures plotting CLvAlpha, CDvAlpha, CL^2vCD, and CLvCD. All
%of these plots have the polyfit regressions that exclude the stall data,
%and the actual values that also exclude the stall data so the plots can be
%seen to more accurate.
figure()
plot(alpha,CL,alpha(1:14),polyval(CLvAlpha,alpha(1:14)))
xlabel('Angle of Attack(\alpha)')
ylabel('CL')
title('\alpha vs. CL')
legend('CFD data','polyfit data')

figure()
plot(alpha_trimmed,CD_trimmed,alpha_trimmed,polyval(CDvAlpha,alpha_trimmed))
xlabel('Angle of Attack(\alpha)')
ylabel('CD')
title('\alpha vs. CD')
legend('CFD data','polyfit data')

figure()
hold on;
plot(CL_trimmed.^2,CD_trimmed,CL_trimmed.^2,polyval(p_CL_2vCD,CL_trimmed.^2))
plot(CL,CD)
title('CL^2 vs CD')
ylabel('CD')
xlabel('CL^2')
legend('CFD data','polyfit data','total data(CLvCD)','location','northwest')

figure()
plot(CL_trimmed,CD_trimmed,CL_trimmed,polyval(p_CLvCD,CL_trimmed))
title('CLvCD')
ylabel('CD')
xlabel('CL')
legend('CFD Data','Polyfit data')


%Below I solve for the k value and CD_0 value by using the equation of
%CD = CD_0 + kCL^2. In this equation we can see that the k is the slope of
%the relationship between CL^2 and CD so creating a linear regression of CD
%vs CL^2 will give me a slope of k. I can also create a quadratic
%regression of CL vs CD and evaluate that at CL=0 to get CD_0.
k = p_CL_2vCD(1)
CD_0 = polyval(p_CLvCD,0) %CD where CL=0

%below are the calculations of gross weight to be used in the max range and
%max endurance calculations, as well as the velocity needed for max range
%and endurance
Wgross = W_0*9.81;
U_E = sqrt((2*Wgross)/(rho*S)*(sqrt(k/(3*CD_0))));
U_R = sqrt((2*Wgross)/(rho*S)*(sqrt(k/(CD_0))));

%Below are the equations used in calculating maximum range and endurance
E_max = ((Rt)^(1-n)) * ((n_tot*Volt*Cap) / ((2/sqrt(rho*S)) * (CD_0^(1/4)) * (2*Wgross*sqrt(k/3))^(3/2)))^n
R_max = (((Rt)^(1-n)) * ((n_tot*Volt*Cap) / ((1/sqrt(rho*S)) * (CD_0^(1/4)) * (2*Wgross*sqrt(k))^(3/2)))^n)*U_R*3.6