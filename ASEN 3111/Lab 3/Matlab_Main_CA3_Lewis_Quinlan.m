%% ASEN 3111 - Computational Assignment 03 - Main
%For this assignment we were tasked with implementing the vortex panel
%method of finding a coeffient of both lift and pressure for different NACA
%airfoiils that all have different amounts of camber as well as thickness.
%In order to do this I had to create functions to evaluate different values
%of x and y that represent the boundary coordinates of the airfoil for the
%coefficient of lift and pressure. I had to make a seperate function to
%actually find these boundary points before I could actually use the vortex
%panel method.

% Author: Quinlan Lewis
% Collaborators: Shawn Stone
% Date: 10/29/20
% Last Modified: 10/29/20

%% Housekeeping
close all; clc; clear;

%% Beginning constants
%these constants were chosen based off of the Kueth & Chow documentation in
%order to test my results from this code with that of the given fortran
%code
c = 1; %[meters] for chord length
N = 1000; %number of panels
ylimits = [-.25 .25];
alpha = deg2rad(0);
V_inf = 60; %[m/s]

%% NACA 0012; Error, and number of nominal panels
m = 0/100;
p = 0/10;
t = 12/100;

%The below function is used throughout the lab in order to calculate the
%boundary points of the different NACA airfoils
[x,y,a_cl_01] = NACA_Airfoils(m,p,t,c,N,'0012',ylimits,"No");

%Vortex panel method used to calculate the coefficient of pressure and
%coefficeint of lift
[~,cp_1] = Vortex_Panel(x,y,V_inf,alpha,(length(x)-1),'0012','No');

%plotting different numbers of panels vs Cp to calculate a nominal number of
%panels to have an error of less that 1% 
variation = 100:50:1000;
cp_1_error = zeros(length(variation),1);
figure(); hold on;
for i = 1:length(variation)
    [x,y,~] = NACA_Airfoils(m,p,t,c,variation(i),'0012',ylimits,"No");
    [~,cp_1_error(i)] = Vortex_Panel(x,y,V_inf,alpha,(length(x)-1),'0012',"Yes");
    cp_1_error(i) = abs((cp_1 - cp_1_error(i))/cp_1) * 100;
end
set(gca,'YDir','reverse')
ylabel('Cp'); xlabel('x/c'); title(sprintf('NACA %s, \\alpha: %d\260, number of panels varying from 198 to 1998','0012',rad2deg(alpha)))
hold off;

%find the index in the Cp array where the percent error is less than 1%
%this value will be the nominal number of panels
index = find(cp_1_error<1,1);
N = variation(index);

%plotting the variation of panels vs the percentage error when using N =
%1000 as the actual value of panels
figure()
plot(variation,cp_1_error);
title('% error of Cp vs Number of panels used'); ylabel('% error'); xlabel('Number of Panels')

%printing nominal number of panels to command window
N_actual = N*2-2;
fprintf('Nominal number of panels is N = %d panels\n',N_actual)

%plotting all given alphas for 0012 airfoil
alphas = deg2rad([-5; 0; 5; 10]);
[x,y,~] = NACA_Airfoils(m,p,t,c,N,'0012',ylimits,"No");
figure(); hold on;
for i = 1:length(alphas)
    [~,~] = Vortex_Panel(x,y,V_inf,alphas(i),(length(x)-1),'0012',"Yes");
end
set(gca,'YDir','reverse')
ylabel('Cp'); xlabel('x/c'); title(sprintf('NACA %s, Varying Alpha, Panels: %d','0012',N_actual))
legend('\alpha: -5','\alpha: 0','\alpha: 5','\alpha: 10')
hold off;

%% Plotting CL slope vs alpha for all Airfoils
%getting cl values for varying alpha between -5 and 15 degrees
alphas = deg2rad(-5:1:20);
cl_1 = zeros(length(alphas),1);
cl_2 = zeros(length(alphas),1);
cl_3 = zeros(length(alphas),1);
cl_4 = zeros(length(alphas),1);

%parameters for 2412 airfoil
m = 2/100;
p = 4/10;
t = 12/100;
[x2,y2,a_cl_02] = NACA_Airfoils(m,p,t,c,N,'2412',ylimits,"No");

%parameters for 4412 airfoil
m = 4/100;
p = 4/10;
t = 12/100;
[x3,y3,a_cl_03] = NACA_Airfoils(m,p,t,c,N,'2412',ylimits,"No");

%parameters for 2424 airfoil
m = 2/100;
p = 4/10;
t = 24/100;
[x4,y4,a_cl_04] = NACA_Airfoils(m,p,t,c,N,'2412',ylimits,"No");

for i = 1:length(alphas)
    [cl_1(i),~] = Vortex_Panel(x,y,V_inf,alphas(i),(length(x)-1),'0012',"Skip");
    [cl_2(i),~] = Vortex_Panel(x2,y2,V_inf,alphas(i),(length(x)-1),'2412',"Skip");
    [cl_3(i),~] = Vortex_Panel(x3,y3,V_inf,alphas(i),(length(x)-1),'4412',"Skip");
    [cl_4(i),~] = Vortex_Panel(x4,y4,V_inf,alphas(i),(length(x)-1),'2424',"Skip");
end

figure(); hold on;
plot(cl_1,alphas); plot(cl_2,alphas); plot(cl_3,alphas); plot(cl_4,alphas)
ylabel('C_l'); xlabel('\alpha'); title('C_l vs. \alpha, for all Airfoils')
yline(0,'--r'); 
legend('NACA 0012','NACA 2412','NACA 4412','NACA 2424','Lift = 0 line','location','Northwest');

%hardcoded values of where the coefficient is closest to 0
a_cl0_vort1 = rad2deg(alphas(6));
a_cl0_vort2 = rad2deg(alphas(4));
a_cl0_vort3 = rad2deg(alphas(2));
a_cl0_vort4 = rad2deg(alphas(4));

%printing the TAT and vortex panel method aoa for zero lift to compare them
fprintf('NACA 0012, TAT alpha_L=0: %f, Vortex Panel alpha_L=0: %f\n',a_cl_01,a_cl0_vort1);
fprintf('NACA 2412, TAT alpha_L=0: %f, Vortex Panel alpha_L=0: %f\n',a_cl_02,a_cl0_vort2);
fprintf('NACA 4412, TAT alpha_L=0: %f, Vortex Panel alpha_L=0: %f\n',a_cl_03,a_cl0_vort3);
fprintf('NACA 2424, TAT alpha_L=0: %f, Vortex Panel alpha_L=0: %f\n',a_cl_04,a_cl0_vort4);
