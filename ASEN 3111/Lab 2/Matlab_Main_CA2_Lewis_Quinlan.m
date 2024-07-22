%% ASEN 3111 - Computational Assignment 01 - Main
%For this assignment, we were tasked with using thin airfoil theory to plot
%the contours of pressure, stream function, and equipotential lines acting
%on a flat plate with a chord length c. We were also tasked with determing
%the relationships between changing the angle of attack, chord length, and
%free stream velocity and how they effect the contours. To do these tasks,
%I used the principle of super position to sum the effects of uniform flow
%and vortex flow acting on a flat airfoil.

% Author: Quinlan Lewis
% Collaborators: Shawn Stone
% Date: 10/8/20
% Last Modified: 10/8/20

%% Housekeeping
clear all; clc; close all;

 %% Declaring Constants
c = 2;
alpha = deg2rad(12);
V_inf = 68;
p_inf = 101.3e3;
rho_inf = 1.225;
N=100;

%% Calculating Phi, Psi, and Pressure for Baseline
%finds psi phi and pressure for N = 100 vortices
[psi_1,phi_1,P_1,V_1] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N,'yes');

%% Calculating Error
%finds psi phi and pressure for N = 10000 vortices to represent values of
%psi phi and pressure that are more accurate that N = 100 panels, this
%value will be used as the actual
N = 1090;
[psi_actual,phi_actual,P_actual,V_actual] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N,'no');

%Preallocating vectors and memmory for the error in Velocity and Pressure
N_error = 100:10:1090;
V_error = zeros(1000/10,1);
P_error = zeros(1000/10,1);
%The below for loop calculates the error in the Pressure and Velocity for
%number of vortices varying from 100 to 1000 and comparing these
%values of pressure and velocity to an actual value of velocity and
%pressure calculated by using number of vortices being 1090.
for i = 1:(1000/10)
    [psi,phi,P,V] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N_error(i),'no');
    
    V_error(i) = mean(abs((V - V_actual)./V_actual)*100,'all');
    P_error(i) = mean(abs((P - P_actual)./P_actual)*100,'all');
end
%Plots the error of pressure and velocity on the same plot
figure()
plot(N_error,V_error,'LineWidth',2)
hold on;
plot(N_error,P_error,'LineWidth',2)
title('Error in Pressure and Velocity')
xlabel('Number of Vortices, N, used')
ylabel('Percentage Error')
legend('Velocity Error','Pressure Error')
hold off;

%% Study On how c effects the function
%reset number of vortices to be 100 for faster calculations
N = 100;
%Changing cord length to vary from 1 to 10 meters
c_vec = 1:3:10;
count = 1;
figure()
for i = 1:length(c_vec)
    c = c_vec(i);
    
    % Define Domain to use contour plot
    xmin=-c-1;
    xmax=c+1;
    ymin=-2;
    ymax=2;
    
    % Define Number of Grid Points for contour plot
    nx=100; % steps in the x direction
    ny=100; % steps in the y direction
    
    % Create mesh over domain using number of grid points specified to be
    % used to plot the functions with contour fucntion
    [x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
    
    %run plot_airfoil_flow function to determine values of the
    %streamfunction, equipotential lines, pressure and velocity
    [psi,phi,P,~] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N_error(i),'no');
    
    %Below the values that were calculated through the plot_airfoil_flow
    %function are added onto a subplot
    subplot(length(c_vec),3,count)
    contour(x,y,psi,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Stream Lines Contour Plot(C = %d)',c_vec(i)))
    count = count + 1;
    
    subplot(length(c_vec),3,count)
    contour(x,y,phi,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Equipotential Lines Contour Plot(C = %d)',c_vec(i)))
    count = count + 1;
    
    subplot(length(c_vec),3,count)
    contourf(x,y,P,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Pressure Contours Plot(C = %d)',c_vec(i)))
    count = count + 1;
end
sgtitle('Study Of how Function is Effected by Chord Length')

%% Study On how alpha effects the function
%Reset chord length to original size
c = 2;
%redefine domain to be used in contour plot
xmin=-c-1;
xmax=c+1;
ymin=-2;
ymax=2;

% Define Number of Grid Points
nx=100; % steps in the x direction
ny=100; % steps in the y direction

% Create mesh over domain using number of grid points specified
[x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
    
%Changing angle of attack
alpha_vec = 8:4:20;
count = 1;
figure()
for i = 1:length(alpha_vec)
    alpha = deg2rad(alpha_vec(i));
    
    %run plot_airfoil_flow function to determine values of the
    %streamfunction, equipotential lines, pressure and velocity with
    %different value for angle of attack
    [psi,phi,P,~] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N_error(i),'no');
    
    %Below the values that were calculated through the plot_airfoil_flow
    %function are added onto a subplot    
    subplot(length(alpha_vec),3,count)
    contour(x,y,psi,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Stream Lines Contour Plot(Alpha = %d)',alpha_vec(i)))
    count = count + 1;
    
    subplot(length(alpha_vec),3,count)
    contour(x,y,phi,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Equipotential Lines Contour Plot(Alpha = %d)',alpha_vec(i)))
    count = count + 1;
    
    subplot(length(alpha_vec),3,count)
    contourf(x,y,P,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Pressure Contours Plot(Alpha = %d)',alpha_vec(i)))
    count = count + 1;
end
sgtitle('Study Of how Function is Effected by Angle of Attack')

%% Study On how alpha effects the function
%resets alpha to be 12 deg
alpha = deg2rad(12);
    
%Changing freestream air velocity acting on the flat plate
Vinf_vec = 40:20:100;
count = 1;
figure()
for i = 1:length(alpha_vec)
    V_inf = Vinf_vec(i);
    
    %run plot_airfoil_flow function to determine values of the
    %streamfunction, equipotential lines, pressure and velocity with
    %different value for freestream air velocity    
    [psi,phi,P,~] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N_error(i),'no');
    
    %Below the values that were calculated through the plot_airfoil_flow
    %function are added onto a subplot     
    subplot(length(Vinf_vec),3,count)
    contour(x,y,psi,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Stream Lines Contour Plot(Vinf = %d)',Vinf_vec(i)))
    count = count + 1;
    
    subplot(length(Vinf_vec),3,count)
    contour(x,y,phi,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Equipotential Lines Contour Plot(Vinf = %d)',Vinf_vec(i)))
    count = count + 1;
    
    subplot(length(Vinf_vec),3,count)
    contourf(x,y,P,100)
    colorbar
    hold on
    plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
    title(sprintf('Pressure Contours Plot(Vinf = %d)',Vinf_vec(i)))
    count = count + 1;
end
sgtitle('Study Of how Function is Effected by Freestream Air Speed')