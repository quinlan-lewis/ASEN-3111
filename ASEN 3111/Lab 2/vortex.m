function [Vort_psi, Vort_phi, u, v] = vortex(x,y,x_v,y_v,G_v)
    %This function takes in the x and y meshgrid, the x and y position of
    %the current vortex, the gamma being used. 
    %The function then calculates the stream function(psi), the
    %equipotential lines(phi), and the u and v components that are
    %generated from the vortices

    r = sqrt((x-x_v).^2 + (y-y_v).^2);
    Vort_psi = (G_v/(2*pi))*log(r);
    
    theta = atan2(-(y-y_v),-(x-x_v));
    Vort_phi = -(G_v/(2*pi))*theta;
    
    %These u and v functions were derived on paper and then checked using
    %an online derivative calculator
    u = G_v*((y-y_v)./((2*pi).*(r.^2)));    
    v = -G_v*((x-x_v)./((2*pi).*(r.^2)));
end

