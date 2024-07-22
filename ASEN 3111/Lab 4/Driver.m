%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Lab #4:   %
% The Vortex Panel Method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1: Implement Vortex Panel Method from Kuethe and Chow %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free-Stream Speed and Angle of Attack

VINF = 10;
ALPHA= 8;

%%%%%%%%%%%%%%%%%%%%%%%
% Kuethe and Chow Data

XB = [1 .933 .750 .500 .250 .067 0 .067 .250 .500 .750 .933 1];
YB = [0 -.005 -.017 -.033 -.042 -.033 0 .045 .076 .072 .044 .013 0];

CP_KC = [0.2630 0.1969 0.2097 0.2667 0.4707 0.9929 -1.8101 -1.5088 -0.9334 -0.5099 -0.1688 0.1674];

%%%%%%%%%%%%%%%%%%%%
% Call Vortex_Panel

[~,CP,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare with Kuethe and Chow Results

fprintf('PROBLEM 1 PRESSURE COEFFICIENT VERIFICATION TEST\nKuethe and Chow (Left) Versus Computed (Right) \n');
disp([CP_KC; CP]');
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2: Use Vortex Panel Method to Study NACA 0012            %
%    Part 1: Determine Suitable Number of Panels                   %
%                                                                  %
%    Approach: Converge Standard Deviation of Pressure Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Airfoil Characteristics

m = 0.00;
p = 0.0;
t = 0.12;
c = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free-Stream Speed and Angle of Attack

VINF = 10;
ALPHA= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Standard Deviation of Pressure Coefficient Using Many Panels

N_half = 1000;
[XB,YB] = NACA_Airfoil(m,p,t,c,N_half);
[CL,CP,XC] = Vortex_Panel(XB,YB,VINF,ALPHA,0);

STD_Actual = std(CP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine # of Panels Required for 5.0% Relative Error

for N_half = 5:1000
    [XB,YB] = NACA_Airfoil(m,p,t,c,N_half);
    [CL,CP,XC] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
    
    STD_Computed = std(CP);
    
    if abs(STD_Actual - STD_Computed)/abs(STD_Actual) < 0.05
        N_half_required = N_half;
        break;
    end
end

fprintf('PROBLEM 2 NUMBER OF PANELS REQUIRED TO ATTAIN STANDARD DEVIATION OF PRESSURE COEFFICIENT WITH 5.0 PERCENT RELATIVE ERROR \n');
disp(2*N_half_required);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2: Use Vortex Panel Method to Study NACA 0012            %
%    Part 2: Plot CL and CP for alpha = -5, 0, 5, 10 Degrees       %
%                                                                  %
%    Approach: Converge Standard Deviation of Pressure Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Construct Boundary Points

[XB,YB] = NACA_Airfoil(m,p,t,c,N_half_required);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CP and CL for Angle of Attack of -5 Degrees

ALPHA = -5;
[CL(1),CP(1,:),XC(1,:)] = Vortex_Panel(XB,YB,VINF,ALPHA,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CP and CL for Angle of Attack of 0 Degrees

ALPHA = 0;
[CL(2),CP(2,:),XC(2,:)] = Vortex_Panel(XB,YB,VINF,ALPHA,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CP and CL for Angle of Attack of 5 Degrees

ALPHA = 5;
[CL(3),CP(3,:),XC(3,:)] = Vortex_Panel(XB,YB,VINF,ALPHA,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CP and CL for Angle of Attack of 10 Degrees

ALPHA = 10;
[CL(4),CP(4,:),XC(4,:)] = Vortex_Panel(XB,YB,VINF,ALPHA,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Sectional Coefficients of Lift

fprintf('PROBLEM 2 NACA 0012 SECTIONAL COEFFICIENTS OF LIFT \n');
fprintf('  AOA of -5 Degrees: %f \n',CL(1));
fprintf('  AOA of 0 Degrees: %f \n',CL(2));
fprintf('  AOA of 5 Degrees: %f \n',CL(3));
fprintf('  AOA of 10 Degrees: %f \n',CL(4));
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Coefficients of Pressure

figure(1)
plot(XC(1,1:N_half_required),-CP(1,1:N_half_required),'--r',XC(1,N_half_required+1:2*N_half_required),-CP(1,N_half_required+1:2*N_half_required),':r','LineWidth',2)
hold on
plot(XC(2,1:N_half_required),-CP(2,1:N_half_required),'--b',XC(2,N_half_required+1:2*N_half_required),-CP(2,N_half_required+1:2*N_half_required),':b','LineWidth',2)
plot(XC(3,1:N_half_required),-CP(3,1:N_half_required),'--g',XC(3,N_half_required+1:2*N_half_required),-CP(3,N_half_required+1:2*N_half_required),':g','LineWidth',2)
plot(XC(4,1:N_half_required),-CP(4,1:N_half_required),'--k',XC(4,N_half_required+1:2*N_half_required),-CP(4,N_half_required+1:2*N_half_required),':k','LineWidth',2)
l = legend('$\alpha = -5^\circ$ Lower','$\alpha = -5^\circ$ Upper','$\alpha = 0^\circ$ Lower','$\alpha = 0^\circ$ Upper','$\alpha = 5^\circ$ Lower','$\alpha = 5^\circ$ Upper','$\alpha = 10^\circ$ Lower','$\alpha = 10^\circ$ Upper');
set(l,'FontName','Times New Roman','FontSize',18,'interpreter','latex')
xlabel('$x/c$','interpreter','latex','FontSize',18)
ylabel('$-C_p$','interpreter','latex','FontSize',18)
set(gca,'FontName','Times New Roman','FontSize',16)
title('Problem 2 NACA 0012 Coefficients of Pressure','interpreter','latex','FontSize',18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3: Plot CL vs AOA for Several Different NACA Airfoils %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attain CL vs AOA for NACA 0012

m = 0.00; p = 0.0; t = 0.12; c = 1;
[XB,YB] = NACA_Airfoil(m,p,t,c,N_half_required);

for ALPHA = -5:1:15
    [CL(1,6+ALPHA),~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
end

ALPHA = -5:1:15;
Fit = polyfit(ALPHA'*pi/180,squeeze(CL(1,:)),1);

VP_Slope(1) = Fit(1);
VP_ZeroLiftAOA(1) = (-Fit(2)/Fit(1))*(180/pi);

TA_Slope(1) = 2*pi;
TA_ZeroLiftAOA(1) = NACA_Airfoil_ZeroLiftAOA(m,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attain CL vs AOA for NACA 2412

m = 0.02; p = 0.4; t = 0.12; c = 1;
[XB,YB] = NACA_Airfoil(m,p,t,c,N_half_required);

for ALPHA = -5:1:15
    [CL(2,6+ALPHA),~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
end

ALPHA = -5:1:15;
Fit = polyfit(ALPHA'*pi/180,squeeze(CL(2,:)),1);

VP_Slope(2) = Fit(1);
VP_ZeroLiftAOA(2) = (-Fit(2)/Fit(1))*(180/pi);

TA_Slope(2) = 2*pi;
TA_ZeroLiftAOA(2) = NACA_Airfoil_ZeroLiftAOA(m,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attain CL vs AOA for NACA 4412

m = 0.04; p = 0.4; t = 0.12; c = 1;
[XB,YB] = NACA_Airfoil(m,p,t,c,N_half_required);

for ALPHA = -5:1:15
    [CL(3,6+ALPHA),~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
end

ALPHA = -5:1:15;
Fit = polyfit(ALPHA'*pi/180,squeeze(CL(3,:)),1);

VP_Slope(3) = Fit(1);
VP_ZeroLiftAOA(3) = (-Fit(2)/Fit(1))*(180/pi);

TA_Slope(3) = 2*pi;
TA_ZeroLiftAOA(3) = NACA_Airfoil_ZeroLiftAOA(m,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attain CL vs AOA for NACA 2424

m = 0.02; p = 0.4; t = 0.24; c = 1;
[XB,YB] = NACA_Airfoil(m,p,t,c,N_half_required);

for ALPHA = -5:1:15
    [CL(4,6+ALPHA),~,~] = Vortex_Panel(XB,YB,VINF,ALPHA,0);
end

ALPHA = -5:1:15;
Fit = polyfit(ALPHA'*pi/180,squeeze(CL(4,:)),1);

VP_Slope(4) = Fit(1);
VP_ZeroLiftAOA(4) = (-Fit(2)/Fit(1))*(180/pi);

TA_Slope(4) = 2*pi;
TA_ZeroLiftAOA(4) = NACA_Airfoil_ZeroLiftAOA(m,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot CL vs AOA for All Airfoils

figure(2)

ALPHA = -5:1:15;
plot(ALPHA,CL(1,:),ALPHA,CL(2,:),ALPHA,CL(3,:),ALPHA,CL(4,:))
l = legend('NACA 0012','NACA 2412','NACA 4412','NACA 2424');
set(l,'FontName','Times New Roman','FontSize',18,'Location','Best')
xlabel('$\alpha$','interpreter','latex','FontSize',18)
ylabel('$c_l$','interpreter','latex','FontSize',18)
set(gca,'FontName','Times New Roman','FontSize',16)
title('Problem 3 NACA Airfoil Sectional Coefficients of Lift Versus AOA','interpreter','latex','FontSize',18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Report Lift Slopes and Zero Lift AOAs

fprintf('PROBLEM 3 NACA Airfoil Lift Slopes and Zero Lift Angles of Attack \n');
fprintf('  NACA 0012\n')
fprintf('     Vortex Panel Lift Slope: %f per Radian \n',VP_Slope(1));
fprintf('     Thin Airfoil Lift Slope: %f per Radian \n',TA_Slope(1));
fprintf('     Vortex Panel Zero Lift AOA: %f Degrees \n',VP_ZeroLiftAOA(1));
fprintf('     Thin Airfoil Zero Lift AOA: %f Degrees \n',TA_ZeroLiftAOA(1));
fprintf('  NACA 2412\n')
fprintf('     Vortex Panel Lift Slope: %f per Radian \n',VP_Slope(2));
fprintf('     Thin Airfoil Lift Slope: %f per Radian \n',TA_Slope(2));
fprintf('     Vortex Panel Zero Lift AOA: %f Degrees \n',VP_ZeroLiftAOA(2));
fprintf('     Thin Airfoil Zero Lift AOA: %f Degrees \n',TA_ZeroLiftAOA(2));
fprintf('  NACA 4412\n')
fprintf('     Vortex Panel Lift Slope: %f per Radian \n',VP_Slope(3));
fprintf('     Thin Airfoil Lift Slope: %f per Radian \n',TA_Slope(3));
fprintf('     Vortex Panel Zero Lift AOA: %f Degrees \n',VP_ZeroLiftAOA(3));
fprintf('     Thin Airfoil Zero Lift AOA: %f Degrees \n',TA_ZeroLiftAOA(3));
fprintf('  NACA 2424\n')
fprintf('     Vortex Panel Lift Slope: %f per Radian \n',VP_Slope(4));
fprintf('     Thin Airfoil Lift Slope: %f per Radian \n',TA_Slope(4));
fprintf('     Vortex Panel Zero Lift AOA: %f Degrees \n',VP_ZeroLiftAOA(4));
fprintf('     Thin Airfoil Zero Lift AOA: %f Degrees \n',TA_ZeroLiftAOA(4));