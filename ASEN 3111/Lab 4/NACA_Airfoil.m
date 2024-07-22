function [XB,YB] = NACA_Airfoil(m,p,t,c,N_half)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                          %
%                                 %
% m = Thickness of Camber Line    %
% p = Location of Max Thickness   %
% t = Thickness                   %
% c = Chord Length                %
% N_half = Half the # of Panels   %
%                                 %
% Output:                         %
%                                 %
% XB = Boundary Points x-location %
% YB = Boundary Points y-location %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(0,c,N_half+1);

y_t = t/0.2*(0.2969*sqrt(x/c)-0.1260*(x/c)-0.3516*(x/c).^2+0.2843*(x/c).^3-0.1036*(x/c).^4);

if m == 0 && p == 0
    for i = 1:N_half+1
       y_c(i) = 0;
       dy_c(i) = 0;
    end
else    
    for i = 1:N_half+1
        if x(i) <= p*c
            y_c(i) = m*x(i)/p^2*(2*p-x(i)/c);
            dy_c(i) = 2*m/p^2*(p-x(i)/c);
        else
            y_c(i) = m*(c-x(i))/(1-p)^2*(1+x(i)/c-2*p);
            dy_c(i) = 2*m/(1-p)^2*(p-x(i)/c);
        end
    end
end

xi = atan(dy_c);

XB_U = x - y_t.*sin(xi);
XB_L = x + y_t.*sin(xi);
YB_U = y_c + y_t.*cos(xi);
YB_L = y_c - y_t.*cos(xi);

XB = [fliplr(XB_L(2:end)) XB_U];
YB = [fliplr(YB_L(2:end)) YB_U];