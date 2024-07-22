syms theta;
eq = tand(17.5) == (2/tand(theta))*((2.5^2*(sind(theta))^2-1)/(2.5^2*(1.4+cosd(2*theta))+2))

vpasolve(eq,theta)