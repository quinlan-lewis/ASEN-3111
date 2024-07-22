function [xout,yout,alpha_L0] =  NACA_Airfoils(m,p,t,c,N,string,ylimits,string2)
%the purpose of this function is to create the boundary points in the x and
%y plane for the given NACA Airfoil entered. This method comes from the
%given lab document where we are given the definitions of the piece wise
%functions to use to find the x and y coordinates of the upper and lower
%boundaries of the airfoil. 

    x = linspace(0,c,N);
    
    xfinal = zeros(length(x),2);
    yfinal = zeros(length(x),2);
    
    for i = 1:length(x)
        y_t = (t/0.2)*c*(0.2969*sqrt(x(i)/c) - 0.1260*(x(i)/c) - 0.3516*(x(i)/c)^2 + 0.2843*(x(i)/c)^3 - 0.1036*(x(i)/c)^4);
        
        if (x(i) >= 0) && (x(i) < p*c)
            y_c = (m*x(i)/p^2)*(2*p - x(i)/c);
            eps = atan(-(m*x(i))*(3*x(i) - 4*c*p)/(c*p));
        elseif (x(i) >= p*c) && (x(i) <= c)
            y_c = m*((c-x(i))/(1-p)^2)*(1 + x(i)/c - 2*p);
            eps = atan(-(2*m*(x(i)-c*p))/(c*(p-1)^2));
        end
                
        x_u = x(i) - y_t*sin(eps);
        x_l = x(i) + y_t*sin(eps);
        
        y_u = y_c + y_t*cos(eps);
        y_l = y_c - y_t*cos(eps);
        
        xfinal(i,:) = [x_u,x_l];
        yfinal(i,:) = [y_u,y_l];
    end
    
    %concatinating the data so the upper and lower data points of the
    %airfoil are in a single array
    xout = cat(1,flip(xfinal(2:length(xfinal),2)),(xfinal(:,1)));
    yout = cat(1,flip(yfinal(2:length(yfinal),2)),(yfinal(:,1)));
    
    %the below if statement is used to plot the airfoils if the user would
    %like to see the geometry of the inputted airfoil
    if string2 == "Yes"
        figure();
        plot(xout,yout)
        title(sprintf('NACA %s',string))
        ylabel('m');
        xlabel('m');
        ylim(ylimits);
    end
    
    %the below if statement is to determine values of the angle of attack
    %where the coefficient of lift is zero
    if p ~= 0
        x = @(theta) (1/pi).*(-1.*(2.*m.*(((c./2).*(1-cos(theta))) - c.*p))./(c.*(p - 1)^2)).*(1 - cos(theta)) + (1/pi).*(-1.*(m.*((c/2).*(1-cos(theta))).*(3.*((c/2).*(1-cos(theta)) - 4.*c.*p)./(c.*p)).*(1 - cos(theta))));
        alpha_L0 = rad2deg(integral(x, 0, pi));
    else
        alpha_L0 = 0;
    end
    
end

