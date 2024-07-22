function [alpha0] = NACA_Airfoil_ZeroLiftAOA(m,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                        %
%                               %
% m = Thickness of Camber Line  %
% p = Location of Max Thickness %
% t = Thickness                 %
%                               %
% Output:                       %
%                               %
% alpha0 = Zero Lift AOA        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = linspace(0,pi,500);
x = (1/2)*(1-cos(theta));

dzdx1 = (2*m/p) * (1 - x/p);
dzdx2 = (-2*m/((1-p)^2)) * (x - p);
dzdx = zeros(size(x));
dzdx(x < p) = dzdx1(x < p);
dzdx(x >= p) = dzdx2(x >= p);

alpha0 = -(1/pi)*trapz(theta,dzdx.*(cos(theta)-1));
alpha0 = alpha0*180/pi;