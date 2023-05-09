%% Turbine influence on velocity profile
% T. Kopperstad 2023

clear; clc; 
close all;

%% Arya Model from Derrick
% Mean Wind Speed Profile in Surface Layer of ABL
% D. Wiberg 2023
% Make sure to have Arya Func and ABLData.xlsx in same folder as this
% function
% Raw data is from Air Force Cambridge Research Lab (AFCRL) experiment in
% Kansas.
% Can play with aero roughness length, z0 (see p 380 of Stull's Boundary Layer
% Meteorology for exact numbers) and displacement distances d to change the
% AFCRL-based profile

% Altitude [m]
z = linspace(0, 120, 120);

% Constants
Beta = 4.7;     % from Businger
k = 0.4;    % Von Karmen constant 

%%%%% Stable BL -  Raw [AFCRL] + Arya %%%%%%

filename = 'ABLData.xlsx';
M = readmatrix(filename,'Sheet','AFCRL','Range','C2:G9');
z0af = 0.0107;             % [From ARCRL]

Z = M(:,1);         % Altitude (m)
U = M(:,3);         % Mean wind speed (m/s)
Theta = M(:,5);      % Potential Temp (K)

l=1; u=8;           % Chooses which 2 data pts to be used in model - the further apart they are the better, generally
[ustar_saf, L_saf, M_saf] = Arya(Z(l), U(l), U(u), Theta(l), Theta(u), k,z0af,z);

%% Wind Turbine variables
rho = 1.225;   % Air density
f = 14.4*2/60; % nominal rpm in hz
hub = 80;      % Hub height [m]
rad = 41;      % Blade radius
A_swept = 5281; % Area swept

a = linspace(0,0.5,11);       % induction ratio
x = linspace(-rad*6,0,100);
x_r = x/rad;
z0 = 0.01;
t0 = (1/f)*log(hub/z0)*rad/hub;
t = linspace(0,100,100);

for i = 1:length(a)
    U_up(i,:) = 1-a(i)*(1+x_r.*(1+x_r.^2).^(-0.5));
    T(i) = 2*pi()*rad^2*M_saf(hub)^2*a(i)*(1-a(i));
    CT(i) = T(i)/(0.5*rho*A_swept*M_saf(hub)^2);
    U_down(i,:) = 0.4*log(t0./t)+CT(i);
end

sym_axis_x = linspace(0,max(M_saf));
sym_axis_y = hub*ones(1,100);

%% Plotting
figure
plot(M_saf, z)
hold on
plot(sym_axis_x,sym_axis_y,"--r")
xlabel('Mean Wind Speed $(m/s)$','Interpreter','latex','Fontsize',14) 
ylabel('Height $(m)$','Interpreter','latex','Fontsize',14)
legend('Arya','Axis of Symmetry','Location','southeast')
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
xlim([0 max(M_saf)])
ylim([0 max(z)])
pbaspect([2,1,1])
%title('Incoming Velocity Profile','Interpreter','latex','Fontsize',36)

figure
plot(x/rad,U_up)
xlabel('Position Upstream of Turbine $(-x/R)$','Interpreter','latex','Fontsize',14)
ylabel('Time averaged velocity $(\bar{u}/\bar{u_\infty})$','Interpreter','latex','Fontsize',14)
%title ('Upstream Velocity Deficit','Interpreter','latex','Fontsize',36)
legend('$a=0$','$a=0.05$','$a=0.1$','$a=0.15$','$a=0.2$','$a=0.25$','$a=0.3$','$a=0.35$','$a=0.4$','$a=0.45$','$a=0.5$','Location','west','orientation','vertical','NumColumns',1,'Interpreter','latex','Fontsize',12)
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
pbaspect([2,1,1])

figure
plot(t,U_down)
xlabel('Transport time $(s)$','Interpreter','latex','Fontsize',14)
ylabel('Time averaged velocity $(\Delta\bar{u}/\bar{u_\infty})$','Interpreter','latex','Fontsize',14)
%title ('Downstream Velocity Deficit','Interpreter','latex','Fontsize',36)
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
%legend('$a=0$','$a=0.05$','$a=0.1$','$a=0.15$','$a=0.2$','$a=0.25$','$a=0.3$','$a=0.35$','$a=0.4$','$a=0.45$','$a=0.5$','Location','eastoutside','orientation','vertical','NumColumns',1,'Interpreter','latex','Fontsize',12)
pbaspect([2,1,1])

figure
plot(a,CT,'o')
hold on
plot(a,CT)
xlabel('Induction factor, $a$','Interpreter','latex','Fontsize',14)
ylabel('Thrust Coefficent $C_T$','Interpreter','latex','Fontsize',14)
%title('Thrust as a function of Induction Factor','Interpreter','latex','Fontsize',36)
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
pbaspect([2,1,1])

%% Roughness Affects on Downstream Velocity

z0 = [0.01 0.05 0.1];
t = linspace(0,100,100);
a = 0.25;
T = 2*pi()*rad^2*M_saf(hub)^2*a*(1-a);
CT = T/(0.5*rho*A_swept*M_saf(hub)^2);
for i = 1:length(z0)
    t0(i) = (1/f)*log(hub/z0(i))*rad/hub;
    U_down_z0(i,:) = 0.4*log(t0(i)./t)+CT;
end

figure
plot(t,U_down_z0)
xlabel('Transport time $(s)$','Interpreter','latex','Fontsize',14)
ylabel('Time averaged velocity $(\Delta\bar{u}/\bar{u_\infty})$','Interpreter','latex','Fontsize',14)
%title ('Downstream Velocity Defecit Roughness variation','Interpreter','latex','Fontsize',36)
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
legend('$z_0=0.01$','$z_0=0.05$','$z_0=0.1$','Location','northeast','orientation','vertical','NumColumns',1,'Interpreter','latex','Fontsize',12)
pbaspect([2,1,1])






