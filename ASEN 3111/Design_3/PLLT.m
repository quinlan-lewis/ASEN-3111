function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,AR)
%% Outputs:
% e = span efficiency factor
% c_L = coefficient of lift
% c_Di = induced coefficient of drag
%% Inputs:
% b = span [ft]
% a0_t = cross-sectional lift slope at tips [per rad]
% a0_r = cross-sectional lift slope at root [per rad]
% c_t = chord at tips [ft]
% c_r = chord at root [ft]
% aero_t = zero-lift angle of attack at the tips [deg]
% aero_r = zero-lift aoa at root [deg]
% geo_t = geo aoa at tip [deg]
% geo_r = geo aoa at root [deg]
% N = number of odd terms to include in series expansion for circulation

%% Calculate constants to be used for CL and CD calculations
%this if statement is used when the aspect ratio is known otherwise the AR
%is passed in as zero as to not change the span.
if AR ~= 0
    b = (AR*(c_r + c_t))/2;
else
    S = (c_r + c_t)*b/2;
    AR = b^2 / S;
end

%% Theta Vector
%these vectors will be used to index into in the PLLT equation so I can
%access the odd terms given as well as the values of theta
i_count = 1:1:(N);
n_odd = 2*i_count - 1;
theta_i = i_count*pi/(2*N);

%% Functions of change in root and tip
%linear functions of span to be used in the PLLT equation
a0_theta = a0_r + (a0_t - a0_r)*cos(theta_i);
c_theta = c_r + (c_t - c_r)*cos(theta_i);
aero_theta = aero_r + (aero_t - aero_r)*cos(theta_i);
geo_theta = geo_r + (geo_t - geo_r)*cos(theta_i);

%% Create system of equations to solve for fourier constants
% preallocate memory for fourier coefficients and assign right hand side of
% equation so I can solve for the unknown fourier coeffs
A = zeros(N,N);
RHS = geo_theta - aero_theta;

%have for loop run to save equations that will be used in fourier coeffs
%calculation
for i = 1:length(theta_i)
    for j = 1:length(theta_i)
        A(i,j) = (4*b/(a0_theta(i)*c_theta(i)))*sin(n_odd(j)*theta_i(i)) + (n_odd(j)*(sin(n_odd(j)*theta_i(i)))/sin(theta_i(i)));
    end
end

%use backslash operator to solve for coeffs
Fourier = A\RHS';

%solve for delta by summation equation
delta = 0;
for i = 2:length(theta_i)
   delta = delta + n_odd(i)*(Fourier(i)/Fourier(1))^2;  
end

%calculate values based on fourier coeffs and delta calculation
e = 1/(1+delta);
c_L = Fourier(1)*pi*AR;
c_Di = (c_L^2)/(pi*AR*e);

end
