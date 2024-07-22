function [c_l,c_p] = Vortex_Panel(x,y,V_inf,alpha,M,string,string2)
%the purpose of this function is to calculate the coefficient of lift and
%coefficient of pressure for the inputted NACA airfoil. The following code
%is based off of the fortran code given to us and originally written by
%Kuethe and Chow.

    N = M+1;
    x_panels = zeros(M,1);
    y_panels = zeros(M,1);
    s = zeros(M,1);
    theta = zeros(M,1);
    sine = zeros(M,1);
    cosine = zeros(M,1);
    RHS = zeros(M,1);
    %the below for loop is made to create data to be used in the actual
    %calculations later in the function
    for i = 1:M
       tmp = i+1;
       x_panels(i) = .5*(x(i)+x(tmp));
       y_panels(i) = .5*(y(i)+y(tmp));
       s(i) = sqrt( (x(tmp)-x(i))^2 + (y(tmp)-y(i))^2);
       theta(i) = atan2( (y(tmp)-y(i)), (x(tmp)-x(i)));
       sine(i) = sin(theta(i));
       cosine(i) = cos(theta(i));
       RHS(i) = sin(theta(i) - alpha);
    end
    
    CN1 = zeros(M,M);
    CN2 = zeros(M,M);
    CT1 = zeros(M,M);
    CT2 = zeros(M,M);
    %the below for loop calculates the matrices CN1, CN2, CT1, and CT2 that
    %will be used to solve a system of equations to get the values for the
    %coefficient of pressure and lift
    for i = 1:M
        for j = 1:M
           if i == j
               CN1(i,j) = -1;
               CN2(i,j) = 1;
               CT1(i,j) = .5*pi;
               CT2(i,j) = .5*pi;
           else
               A = -(x_panels(i)-x(j))*cosine(j) - (y_panels(i)-y(j))*sine(j);
               B = (x_panels(i)-x(j))^2 + (y_panels(i)-y(j))^2;
               C = sin(theta(i)-theta(j));
               D = cos(theta(i)-theta(j));
               E = (x_panels(i)-x(j))*sine(j) - (y_panels(i)-y(j))*cosine(j);
               F = log(1 + s(j)*(s(j)+2*A)/B);
               G = atan2(E*s(j), (B+A*s(j)));
               P = (x_panels(i)-x(j))*sin(theta(i) - 2*theta(j)) + (y_panels(i)-y(j))*cos(theta(i) - 2*theta(j));
               Q = (x_panels(i)-x(j))*cos(theta(i) - 2*theta(j)) - (y_panels(i)-y(j))*sin(theta(i) - 2*theta(j));
               CN2(i,j) = D + .5*Q*F/s(j) - (A*C+D*E)*G/s(j);
               CN1(i,j) = .5*D*F + C*G - CN2(i,j);
               CT2(i,j) = C + .5*P*F/s(j) + (A*D-C*E)*G/s(j);
               CT1(i,j) = .5*C*F - D*G - CT2(i,j);
           end
        end 
    end
    
    AN = zeros(N,N);
    AT = zeros(M,N);
    %The below for loop calculates the matrices AN, and AT to be used in
    %solving the given created system of equations to solve for a value of
    %gamma prime
    for i = 1:M
        AN(i,1) = CN1(i,1);
        AN(i,N) = CN2(i,M);
        AT(i,1) = CT1(i,1);
        AT(i,N) = CT2(i,M);
        for j = 2:M
            AN(i,j) = CN1(i,j) + CN2(i,j-1);
            AT(i,j) = CT1(i,j) + CT2(i,j-1);
        end
    end
    
    AN(N,1) = 1;
    AN(N,N) = 1;
    for j = 2:M
        AN(N,j) = 0;
    end
    RHS(N) = 0;
    
    %this value of gamma prime represents the strength of circulation of
    %the vortices
    gamma_prime = AN\RHS;
    
    %this last for loop calculates the coefficient of pressure at different
    %points along the airfoil
    V = zeros(M,1);
    Cp = zeros(M,1);
    for i = 1:M
        V(i) = cos( theta(i) - alpha);
        for j = 1:N
            V(i) = V(i) + AT(i,j)*gamma_prime(j);
            Cp(i) = 1 - (V(i))^2;
        end
    end
    
    %calculating the cl and average cp which will be used for error
    %analysis
    Gamma = sum(V.*s);
    c_l = 2*Gamma/(max(x) - min(x));
    c_p = mean(Cp);
    
    %the below if statement is used for plotting the coefficient of
    %pressure vs x/c
    if string2  == "Yes"
        plot(x_panels/(max(x) - min(x)),Cp)
    elseif string2 == "Skip"
        return
    else
        figure()
        plot(x_panels/(max(x) - min(x)),Cp)
        set(gca,'YDir','reverse')
        ylabel('Cp'); xlabel('x/c'); title(sprintf('NACA %s, \\alpha: %d\260, Panels: %d',string,rad2deg(alpha),M))
    end

    
end

