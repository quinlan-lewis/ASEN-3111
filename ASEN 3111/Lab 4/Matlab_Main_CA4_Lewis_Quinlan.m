%% ASEN 3111 - Computational Assignment 04 - Main
%For this assignment we were tasked with implementing Prandlt Lifting Line
%Theory to calculate the coefficients of lift, drag, and the span
%efficiency factor. To do this the equation given to us in lecture was
%implemented in the PLLT function. This function calculates all of the
%Fourier coefficients needed for calculating lift and drag coefficients as
%well as e. This function takes in values for span, geometric and
%aerodynamic twist, as well as the number of odd terms that is wanted for
%the calculation. The function is called, then using the returned values for
%the coefficients, the lift and induced drag is calculated.

%NOTE: To complete this code vortex panel method was needed to get the
%angle of attack at lift = 0, as well as the a0 values of the wings. To
%find these values I used the solution code given to us for CA Lab 3
%because I was not confident in my own code from CA lab 3. The functions
%used for those values are "Vortex_Panel.m" and "NACA_Airfoil.m".

% Author: Quinlan Lewis
% Collaborators: Professor Evans(Vortex Panel Method)
% Date: 11/19/20
% Last Modified: 11/19/20

clear; clc; close all;

N_vort = 250;
VINF = 150/2.237;

%% Obtaining a0 and alpha_L=0 for 0012
%CA Solution code is used to get values for VP_Slope(1) and
%VP_ZeroLiftAOA(1) for the tip values

m = 0.00; p = 0.0; t = 0.12; c = 1;
[XB,YB] = NACA_Airfoil(m,p,t,c,N_vort);

for ALPHA = -5:1:15
    [CL(1,6+ALPHA),~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
end

ALPHA = -5:1:15;
Fit = polyfit(ALPHA'*pi/180,squeeze(CL(1,:)),1);

VP_Slope(1) = Fit(1);
VP_ZeroLiftAOA(1) = (-Fit(2)/Fit(1))*(180/pi);

%% Obtaining a0 and alpha_L=0 for 2412
%CA Solution code is used to get values for VP_Slope(2) and
%VP_ZeroLiftAOA(2) for the root values 
m = 0.02; p = 0.4; t = 0.12; c = 1;
[XB,YB] = NACA_Airfoil(m,p,t,c,N_vort);

for ALPHA = -5:1:15
    [CL(2,6+ALPHA),~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
end

ALPHA = -5:1:15;
Fit = polyfit(ALPHA'*pi/180,squeeze(CL(2,:)),1);

VP_Slope(2) = Fit(1);
VP_ZeroLiftAOA(2) = (-Fit(2)/Fit(1))*(180/pi);

%% Constants to be used in PLLT code
%The below constants are given in the lab document and describe the twist
%in the wing as well as the shape of the wing as a whole. Number of odd
%terms is chosen to be 100 to get a close value as to what the actual lift
%and induced drag coeffs should be

b = 100; %[ft]
c_r = 15; %[ft] root chord
c_t = 5; %[ft] tip chord
a0_r = VP_Slope(2); %[per rad]
a0_t = VP_Slope(1); %[per rad]
aero_r = deg2rad(VP_ZeroLiftAOA(2)); %[rad] zero lift aoa
aero_t = deg2rad(VP_ZeroLiftAOA(1)); %[rad] zero lift aoa
geo_r = deg2rad(5); %[rad] geo aoa
geo_t = deg2rad(0); %[rad] geo aoa
N = 100; % number of odd terms

%PLLT function is ran to get e and the corresponding coeffs
[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,0);

%% Lift and Induced Drag acting on A/C
%these values are used in calculating the actual lift and induced drag
%acting on the aircraft

v = 150/2.237; %[m/s]
rho = 1.225;

L_actual = (1/2)*rho*(v^2)*c_L;
Di_actual = (1/2)*rho*(v^2)*c_Di;
fprintf('Lift produced in Newtons: %f\n',L_actual)
fprintf('Induced Drag produced in Newtons: %f\n',Di_actual)

%% Problem 2
%For this problem a series of while loops were created to run until the
%desired error value was met then that corresponding value for the number
%of odd terms was then printed to the command window.

% find number of odd terms so that error is less than 5%.
L_error = 10;
Di_error = 10;
N = 1;
while (L_error > 5) && (Di_error > 5)
    [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,0);
    L_tmp = (1/2)*rho*(v^2)*c_L;
    Di_tmp = (1/2)*rho*(v^2)*c_Di;
    
    L_error = abs((L_actual - L_tmp)/L_actual) * 100;
    Di_error = abs((Di_actual - Di_tmp)/Di_actual) * 100;
    
    N = N + 1;
end
fprintf('Number of odd terms needed for 5 percent relative error: %d\n',N)

%error less than 1%
L_error = 10;
Di_error = 10;
N = 1;
while (L_error > 1) && (Di_error > 1)
    [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,0);
    L_tmp = (1/2)*rho*(v^2)*c_L;
    Di_tmp = (1/2)*rho*(v^2)*c_Di;
    
    L_error = abs((L_actual - L_tmp)/L_actual) * 100;
    Di_error = abs((Di_actual - Di_tmp)/Di_actual) * 100;
    
    N = N + 1;
end
fprintf('Number of odd terms needed for 1 percent relative error: %d\n',N)

%error less than 1/10%
L_error = 10;
Di_error = 10;
N = 1;
while (L_error > 1/10) && (Di_error > 1/10)
    [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,0);
    L_tmp = (1/2)*rho*(v^2)*c_L;
    Di_tmp = (1/2)*rho*(v^2)*c_Di;
    
    L_error = abs((L_actual - L_tmp)/L_actual) * 100;
    Di_error = abs((Di_actual - Di_tmp)/Di_actual) * 100;
    
    N = N + 1;
end
fprintf('Number of odd terms needed for 1/10 percent relative error: %d\n',N)

%% Problem 3
%For this problem we were tasked with plotting the change in span
%efficiency factor(e) with respect to the change in taper ratio. Here the
%taper ratio is defined as c_t/c_r. To answer this problem, I created an
%array of taper ratio values and set my c_t to an arbitrary number. This
%way I can calculate the corresponding c_r value based on the current taper
%ratio and pass those values into the PLLT function.

b = 100; %[ft]
c_r = 0; %[ft] root chord
c_t = 1; %[ft] tip chord
a0_r = 2*pi; %[per rad]
a0_t = 2*pi; %[per rad]
aero_r = 0; %[deg] zero lift aoa
aero_t = 0; %[deg] zero lift aoa
geo_r = deg2rad(5); %[deg] geo aoa
geo_t = deg2rad(5); %[deg] geo aoa
N = 20; % number of odd terms
AR = [4;6;8;10];
ct_o_cr = 0.01:.01:1;

%the below for loop is ran for all of the given values of the aspect ratio
for i = 1:length(AR)
   for j = 1:length(ct_o_cr)
       c_r = c_t/(ct_o_cr(j));
       [e(i,j),~,~] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,AR(i));
   end
end

%all of the vectors of e are plotted on the same graph to show their
%relationship
figure(); hold on;
for i = 1:length(AR)
    plot(ct_o_cr,e(i,:))
end
xlabel('c_t/c_r')
ylabel('e')
title('Span Efficieny factor vs Taper Ratio')
legend('AR = 4','AR = 6','AR = 8','AR = 10','Location','SouthEast')