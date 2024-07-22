%% ASEN 3111 - Computational Assignment 01 - Main
%For this assignment, the task was to use both composite trapezoidal rule
%and zomposite simpsons rule to determine values of coefficient of lift and
%drag as well as values for lift and drag per unit span. The other task was
%to determine how many panels or interation points were required in order
%to reach a certain value for the relative error when compared to the
%actual value. The purpose of this assignment was to understand different
%ways to integrate different functions that sometimes cannot be done
%analytically but instead need to be done numerically.

% Author: Quinlan Lewis
% Collaborators: Shawn Stone
% Date: 9/17/20
% Last Modified: 9/17/20

%% Housekeeping
clear; clc; close all;

%% Problem 1

%% Analytical Solution for Lift and Drag Coefficient
%Solving for coefficient of lift and drag numerically by using symbolic
%functions in matlab
syms R V_inf theta

Gamma = 2*pi*R*V_inf;

C_p = 1 - (4*(sin(theta))^2 + 2*Gamma*sin(theta)/(pi*R*V_inf) + (Gamma/(2*pi*R*V_inf))^2);

cl_calc = double((-1/2)*int(C_p*sin(theta), theta, 0, 2*pi));
cd_calc = double((-1/2)*int(C_p*cos(theta), theta, 0, 2*pi));
fprintf('Lift coefficient, c_l: %.4f\nDrag coefficient, c_d: %.4f\n', (cl_calc), (cd_calc));

%% Calculation of Lift and Drag coefficient using Composite Trapeziodaal Rule
%Using Trapeziodal rule to solve for coeficeint of lift and drag based on
%a variation of number of panels used
panels = 10;
start = 0;
finish = 2*pi;
CL_trap = zeros(1,panels);
CD_trap = zeros(1,panels);
%outside for loop changes the number of panels used in the summation
%portion of the calculation, this also controls the step size that is to be
%used for the given number of panels
for N = 1:panels
    cl = 0;
    cd = 0;
    h = (finish - start)/N;
    th = 0:h:2*pi;
    %embedded for loop calculates sum of however many panels there are for
    %this current iteration
    for i = 1:N
        cl = (cl + ((-1/2)*subs(C_p,theta,th(i+1))*sin(th(i+1)) + (-1/2)*subs(C_p,theta,th(i))*sin(th(i)))/2);
        cd = (cd + ((-1/2)*subs(C_p,theta,th(i+1))*cos(th(i+1)) + (-1/2)*subs(C_p,theta,th(i))*cos(th(i)))/2);
    end
    %calculated cl and cd are then inputted to this vector along with their
    %corresponding number of panels
    CL_trap(N) = double(h*cl);
    CD_trap(N) = double(h*cd);
end

%After calculating the coefficients of lift and drag I plotted them versus
%how many panels are used in finding those values
figure
subplot(2,1,1)
    plot(1:panels,CL_trap);
    xlabel('Number of Panels, N')
    ylabel('Coefficient of Lift, CL')
subplot(2,1,2)
    plot(1:panels,CD_trap);
    xlabel('Number of Panels, N')
    ylabel('Coefficient of Drag, CD')
sgtitle('Composite Trapeziodal Rule for calculating CL and CD')

%% Calculation of Lift and Drag coefficient using Composite Simpson's Rule
%Using Simpson's Rule to solve for coefficient of lift and drag based on a
%variation of number of panels used
CL_simp = zeros(1,panels);
CD_simp = zeros(1,panels);
%outside for loop changes the number of panels used in the summation
%portion of the calculation, this also controls the step size that is to be
%used for the given number of panels
for N = 1:panels 
    cl = 0;
    cd = 0;
    h = (finish - start)/(2*N);
    th = 0:h:2*pi;
    %embedded for loop calculates sum of however many panels there are for
    %this current iteration
    for i = 1:N
        cl = cl + ((-1/2)*subs(C_p,theta,th(2*i-1))*sin(th(2*i-1)) + 4*(-1/2)*subs(C_p,theta,th(2*i))*sin(th(2*i)) + (-1/2)*subs(C_p,theta,th(2*i+1))*sin(th(2*i+1)));
        cd = cd + ((-1/2)*subs(C_p,theta,th(2*i-1))*cos(th(2*i-1)) + 4*(-1/2)*subs(C_p,theta,th(2*i))*cos(th(2*i)) + (-1/2)*subs(C_p,theta,th(2*i+1))*cos(th(2*i+1)));
    end
    %calculated cl and cd are then inputted to this vector along with their
    %corresponding number of panels
    CL_simp(N) = double(h/3*cl);
    CD_simp(N) = double(h/3*cd);
end

%After calculating the coefficients of lift and drag I plotted them versus
%how many panels are used in finding those values
figure
subplot(2,1,1)
    plot(1:panels,CL_simp);
    xlabel('Number of Panels, N')
    ylabel('Coefficient of Lift, CL')
subplot(2,1,2)
    plot(1:panels,CD_simp);
    xlabel('Number of Panels, N')
    ylabel('Coefficient of Drag, CD')
sgtitle('Simpsons Rule for calculating CL and CD')

%% Lift Coefficient relative error calculation
%Calulating how many panels are required in order for trapeziodal rule to
%be accurate to within 1/10 percent
eT = 1;
i = 0;
while eT > (1/10)/100
    i = i+1;
    eT = abs((cl_calc - CL_trap(i))/cl_calc);
end
fprintf('Number of panels required for trapezoidal rule with 1/10%% relative error: %d\n', i);

%Calulating how many panels are required in order for simpsons rule to
%be accurate to within 1/10 percent
eT = 1;
i = 0;
while eT > (1/10)/100
    i = i+1;
    eT = abs((cl_calc - CL_simp(i))/cl_calc);
end
fprintf('Number of panels required for Simpsons rule with 1/10%% relative error: %d\n', i);

%% Problem 2

% assigning all of the constant and given values
load Cp.mat;
c = 2; %chord length [m]
alpha = 9; %angle of attack [deg]
V_inf = 30; %freestream airspeed [m/s]
rho_inf = 1.225; %air density [kg/m^3]
p_inf = 101.3*(10^3); %ambient pressure [Pa]
q_inf = (1/2)*rho_inf*(V_inf)^2;

%% Calculating actual value for Lift and Drag using Composite Trapezoidal rule
%The below for loop finds a more accurate value for the Lift and Drag forces
%acting on the airfoil by using composite trapeziodal rule for 40,000 panels,
%the value found fomr this for loop will be used as the actual value for
%the Lift and Drag per unit span. In this calculation shear stress was left
%out of the integration process because we are assuming an ideal flow.
N_calc = 40000;
N_force = 0;
A_force = 0;
h = (c)/N_calc;
th = 0:h:c;
for i = 1:N_calc
   %calulating the y values from the given function yt by passing the
   %current and next value for x into a seperate function to calculate the
   %corresponding y values
   y = y_t(th(i));
   y1 = y_t(th(i+1));

   %calculates pressure on the upper and lower surface of the airfoil for
   %the current x value and the next x value
   p_u = q_inf*(fnval(Cp_upper, th(i)/c)) + p_inf;
   p_l = q_inf*(fnval(Cp_lower, th(i)/c)) + p_inf;
   p_u1 = q_inf*(fnval(Cp_upper, th(i+1)/c)) + p_inf;
   p_l1 = q_inf*(fnval(Cp_lower, th(i+1)/c)) + p_inf;

   %calulates and sums the normal forces acting on the upper and lower
   %surfaces
   N_force_upper = (-1)*((th(i+1))-(th(i)))*(p_u1 + p_u)/2;
   N_force_lower = ((th(i+1))-(th(i)))*(p_l1 + p_l)/2;
   N_force = N_force + (N_force_upper + N_force_lower);
   
   %calulates and sums the Axial forces acting on the upper and lower
   %surfaces
   A_force_upper = -1*(-1)*(y1 - y)*((p_u1 + p_u)/2);
   A_force_lower = (y1 - y)*((p_l1 + p_l)/2);
   A_force = A_force + (A_force_upper + A_force_lower);
end
Normal_calc = N_force;
Axial_calc = A_force;
L_calc = N_force*cosd(alpha) - A_force*sind(alpha);
D_calc = N_force*sind(alpha) + A_force*cosd(alpha);
fprintf('Lift per unit span of NACA 0012 airfoil: %.2f Newtons/meter\n',L_calc);
fprintf('Drag per unit span of NACA 0012 airfoil: %.2f Newtons/meter\n',D_calc);

%% Calculating the different relative errors for Lift
%For this calculation I orginally ran one whole for loop with 3 different
%if statements, one for each of the necessary relative errors. It ran for
%1:1000 panels, this took about 4 minutes to complete so once that original
%for loop ran I found the correct number of panels needed for each of the
%errors. Using these number of panels, I created 2 seperate for loops, one
%for the 5% and 1% relative error cases and the other for the (1/10)%
%relative error case.

%The for loop below finds values for the Lift and Drag acting on the
%airfoil by composite trapeziodal rule. It runs for an amount of panels
%from 1:200 or any other value of maximum panels the user desires.
panels = 200;
Normal = zeros(1,panels);
Axial = zeros(1,panels);
L = zeros(1,panels);
D = zeros(1,panels);
e1 = 1;
e2 = 1;
for N = 1:panels
    %These declarations reset the values for normal and axial force back to
    %0 for the next number of panels to run and reassign a new step size,
    %h, based on the new number of panels.
    N_force = 0;
    A_force = 0;
    h = (c)/N;
    th = 0:h:c;
    %The below for loop sums all the axial and normal forces of the upper
    %and lower surfaces based on the current number of panels being used.
    for i = 1:N
       %calulating the y values from the given function yt by passing the
       %current and next value for x into a seperate function to calculate the
       %corresponding y values
       y = y_t(th(i));
       y1 = y_t(th(i+1));
       
       %calculates pressure on the upper and lower surface of the airfoil for
       %the current x value and the next x value
       p_u = q_inf*(fnval(Cp_upper, th(i)/c)) + p_inf;
       p_l = q_inf*(fnval(Cp_lower, th(i)/c)) + p_inf;
       p_u1 = q_inf*(fnval(Cp_upper, th(i+1)/c)) + p_inf;
       p_l1 = q_inf*(fnval(Cp_lower, th(i+1)/c)) + p_inf;
       
       %calulates and sums the normal forces acting on the upper and lower
       %surfaces
       N_force_upper = (-1)*((th(i+1))-(th(i)))*(p_u1 + p_u)/2;
       N_force_lower = ((th(i+1))-(th(i)))*(p_l1 + p_l)/2;
       N_force = N_force + (N_force_upper + N_force_lower);
       
       %calulates and sums the Axial forces acting on the upper and lower
       %surfaces       
       A_force_upper = -1*(-1)*(y1 - y)*((p_u1 + p_u)/2);
       A_force_lower = (y1 - y)*((p_l1 + p_l)/2);
       A_force = A_force + (A_force_upper + A_force_lower);
    end
    Normal(N) = N_force;
    Axial(N) = A_force;
    %Calculates the Lift and Drag forces based off of the angle of attack
    %alpha, and the normal and axial forces
    L(N) = N_force*cosd(alpha) - A_force*sind(alpha);
    D(N) = N_force*sind(alpha) + A_force*cosd(alpha);
    
    %This portion determines if the relative error has reached the desired
    %threshhold yet, if it hasn't then nothing will be printed to the
    %command window, if it has then that number of integration points will
    %be printed to the command window.
    error = abs(L_calc - L(N))/L_calc;
    if (error < (5/100) && e1 == 1)
       e1 = 0;
       fprintf('Integration points required for lift solution with 5%% relative error: %d\n',2*(N+1))
    end
    if (error < (1/100) && e2 == 1)
       e2 = 0;
       fprintf('Integration points required for lift solution with 1%% relative error: %d\n',2*(N+1))
    end
end

%this for loop is for calculating the number of integration points required
%for the relative error to be within (1/10)% of the actual value.
panels = 975;
e3 = 1;
for N = 950:panels
    %These declarations reset the values for normal and axial force back to
    %0 for the next number of panels to run and reassign a new step size,
    %h, based on the new number of panels.
    N_force = 0;
    A_force = 0;
    h = (c)/N;
    th = 0:h:c;
    %The below for loop sums all the axial and normal forces of the upper
    %and lower surfaces based on the current number of panels being used.
    for i = 1:N
       %calulating the y values from the given function yt by passing the
       %current and next value for x into a seperate function to calculate the
       %corresponding y values
       y = y_t(th(i));
       y1 = y_t(th(i+1));
       
       %calculates pressure on the upper and lower surface of the airfoil for
       %the current x value and the next x value
       p_u = q_inf*(fnval(Cp_upper, th(i)/c)) + p_inf;
       p_l = q_inf*(fnval(Cp_lower, th(i)/c)) + p_inf;
       p_u1 = q_inf*(fnval(Cp_upper, th(i+1)/c)) + p_inf;
       p_l1 = q_inf*(fnval(Cp_lower, th(i+1)/c)) + p_inf;
       
       %calulates and sums the normal forces acting on the upper and lower
       %surfaces
       N_force_upper = (-1)*((th(i+1))-(th(i)))*(p_u1 + p_u)/2;
       N_force_lower = ((th(i+1))-(th(i)))*(p_l1 + p_l)/2;
       N_force = N_force + (N_force_upper + N_force_lower);
       
       %calulates and sums the Axial forces acting on the upper and lower
       %surfaces       
       A_force_upper = -1*(-1)*(y1 - y)*((p_u1 + p_u)/2);
       A_force_lower = (y1 - y)*((p_l1 + p_l)/2);
       A_force = A_force + (A_force_upper + A_force_lower);
    end
    Normal(N) = N_force;
    Axial(N) = A_force;
    %Calculates the Lift and Drag forces based off of the angle of attack
    %alpha, and the normal and axial forces
    L(N) = N_force*cosd(alpha) - A_force*sind(alpha);
    D(N) = N_force*sind(alpha) + A_force*cosd(alpha);
    
    %This portion determines if the relative error has reached the desired
    %threshhold yet, if it hasn't then nothing will be printed to the
    %command window, if it has then that number of integration points will
    %be printed to the command window.
    error = abs(L_calc - L(N))/L_calc;
    if (error < ((1/10)/100) && e3 == 1)
       e3 = 0;
       fprintf('Integration points required for lift solution with (1/10)%% relative error: %d\n',2*(N+1))
    end
end

%% Exterior functions

%Functions runs for calculating values of y based off of x/c
function y = y_t(x)
%function takes in value for x and plugs into given equation for y
%Author: Quinlan Lewis
    c = 2;
    t = 12/100;
    y = ((t)/.2)*c*(0.2969*sqrt(x/c) - .1260*(x/c) - .3516*(x/c)^2 + .2843*(x/c)^3 - .1036*(x/c)^4);
end

