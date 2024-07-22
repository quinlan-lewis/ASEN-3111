function [StreamFunction,Equipotential,P,V] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N,to_plot)
    %xy values are in cartesian coordinates
    %x = rcos(theta)
    %y = rsin(theta)
    %r = norm([x,y])
    %theta = atan(y/x)
    
    % Define Domain for meshgrid
    xmin=-c-1;
    xmax=c+1;
    ymin=-2;
    ymax=2;

    % Define Number of Grid Points to be on meshgrid
    nx=100; % steps in the x direction
    ny=100; % steps in the y direction
    
    % Create mesh over domain using number of grid points specified
    [x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
    
    %% StreamFunction, Equipotential Lines, and Pressure on Plate
    
    %initializes the streamfunction as only the uniform flow in the field
    uni_psi = V_inf*(cos(alpha).*y - sin(alpha).*x);
    StreamFunction = uni_psi;
    
    %Initializes the equipotential lines to be only effected from uniform
    %flow
    uni_phi = V_inf*(cos(alpha).*x + sin(alpha).*y);
    Equipotential = uni_phi;
    
    %Initializes u and v from uniform flow to add component wise to vortex
    %velocities 
    u = V_inf*cos(alpha);
    v = V_inf*sin(alpha);
    
    %change in distance between each vortex
    delx = c/N;
    
    %initializes the size of the chord in the x and y axis
    ychord = zeros(1,N);
    xchord = linspace(delx/2,(c-delx/2),N);
    
    %Calculates Gamma to be used in finding the effect of the vortex acting
    %on the flow
    gamma = 2*alpha*V_inf.*sqrt(((1-(xchord./c))./(xchord./c)));
    Gamma = gamma.*delx;
    
    %The below for loop adds all of the effects of each vortex into the
    %flow and sums each individual vortex effect to the total flow of the
    %stream function
    for i = 1:N
        [vort_psi,vort_phi,uvort,vvort] = vortex(x,y,xchord(i),ychord(i),Gamma(i));
        StreamFunction = StreamFunction + vort_psi;
        Equipotential = Equipotential + vort_phi;
        u = u + uvort;
        v = v + vvort;        
    end

    %After the velocity components are calculated, they are then summed to
    %find V, which is then used in calculating the pressure acting on the
    %plate
    V = sqrt(u.^2 + v.^2);
    q_inf = 0.5*rho_inf*V_inf^2;
    Cp = 1 - (V./V_inf).^2;
    P = p_inf + q_inf.*Cp;
    
    %% Plotting contours
    %The below plotting only occurs if a conditional statement is added at
    %the function call of this function, that way less plots are made than
    %needed
    if to_plot == "yes"
        figure()
        %plotting stream function contours
        subplot(3,1,1)
        contour(x,y,StreamFunction,100)
        colorbar;
        hold on
        plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
        title(sprintf('Stream Lines Contour Plot(N = %d)',N))

        %plotting equipotential lines contour
        subplot(3,1,2)
        contour(x,y,Equipotential,100)
        colorbar;
        hold on
        plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
        title(sprintf('Equipotential Lines Contour Plot(N = %d)',N))

        %plotting pressure contours
        subplot(3,1,3)
        contourf(x,y,P,100)
        colorbar;
        hold on
        plot([0:c],zeros(size([0:c])),'k','LineWidth',2)
        title(sprintf('Pressure Contours Plot(N = %d)',N))
        
        sgtitle(sprintf('Plots using N = %d',N));
    end
    
end

