    %% Lifting Line Function
clear; clc; 
close all;

% z0 = 0.025;    % Aerodynamic Roughness length [Stull]

% Wind Turbine variables
hub = 80;      % Hub height [m]
rad = 41;      % Blade radius


Beta = 4.7;
k = 0.4;    % Von Karmen constant

z = logspace(-2, 2.5);


% Input Arguments   

filename = 'ABLData.xlsx';

%%%%% Stable BL -  Raw + Arya %%%%%%
% [AFCRL]
M = readmatrix(filename,'Sheet','AFCRL','Range','C2:G9');
z0af = 0.0107;

Zaf = M(:,1);         % Altitude (m)
Uaf = M(:,3);         % Mean wind speed (m/s)
Thetaaf = M(:,5);      % Potential Temp (K)

l=1; u=8;
[ustar_saf, L_saf, M_saf] = Arya(Zaf(l), Uaf(l), Uaf(u), Thetaaf(l), Thetaaf(u), k,z0af,z);


% [Wangara]
z0w = 0.0012;
M = readmatrix(filename,'Sheet','W','Range','C2:H5');

Zw = M(:,1);         % Altitude (m)
Uw = M(:,3);         % Mean wind speed (m/s)
Tw = M(:,5);         % Temp (K)
Pw = M(:,6);
Thetaw = Tw.*(Pw/1000).^0.286;

l=1; u=4;
[ustar_sw, L_sw, M_sw] = Arya(Zw(l), Uw(l), Uw(u), Thetaw(l), Thetaw(u), k,z0w,z);

% ustar_n = k*Uw(l)/log(Zw(l)/z0w);      %[Stull p.383]
% 
% M_n = ustar_n.* (1/k) .* log(z./z0w);

% Simulated
Zi = [1 180];
Ui = [2 14];
Thetai = [293.2900 298.4000];

l=1; u=2;

% Z0 Study
z0i = [.002 .005 .01 .02 .05 .1 .2];

figure
for i = 1:length(z0i)
    z0 = z0i(i);
    [ustar_sz, L_sz, M_sz] = Arya(Zi(l), Ui(l), Ui(u), Thetai(l), Thetai(u), k,z0,z);
    
    semilogy(M_sz, z)
    xstr(i) = strcat("$z_0$ =", string(z0));
    hold on
end

xlabel('Mean Wind Speed [m/s]'); ylabel('Height [m]');
set(gca,'FontSize',14)
set(0,'defaultTextInterpreter','latex');
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
% title('Simulated $z_0$ - Stable')
legend(xstr)

% D Study
z0 = 0.01;
D = [0.5 1 2 5 10];
xstr = strings(length(D),1);


figure
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
for i = 1:length(D)
    d= D(i);
    [ustar_sd, L_sd, M_sd] = Arya(Zi(l), Ui(l), Ui(u), Thetai(l), Thetai(u), k,z0,z-d);
    
    h = d*ones(1,length(z));
    x = linspace(0,length(z)-1,length(z));

    T = semilogy(M_sd, z,x,h, "--");
    T(1).Color = colors{i};
    T(2).Color = colors{i};

    xstr(2*i-1) = strcat("M, d=", string(d));
    xstr(2*i) = strcat("d =", string(d));
    hold on
end
% figure
%     semilogy(M_sd, z,"red", x, h, "--")
%     xstr(i) = strcat("d =", string(d));


xlabel('Mean Wind Speed [m/s]','interpreter','latex');
ylabel('Height [m]','interpreter','latex');
set(gca,'FontSize',14)
% title('Simulated d- Stable','interpreter','latex')
legend(xstr)
xlim([0, max(M_sd)/2])

%%%%% Neutral BL  %%%%%%%
% ustar_n = k*Uaf(8)/log(Zaf(8)/z0af);      %[Stull p.383]
% 
% M_n = ustar_n.* (1/k) .* log(z./z0af);




% Wind Turbine Hub Height
h = hub*ones(1,length(z));
x = linspace(0,length(z)-1,length(z));





% figure
% semilogy(M_n, z, M_s, z, x,h, '--')
% xlabel('Mean Wind Speed [m/s]'); ylabel('Height [m]');
% legend('Neutral','Stable Surface Layer', 'Unstable Surface Layer','Location','southeast');
% xlim([0,10]);
% set(gca,'FontSize',14)

figure
semilogy(M_saf, z, Uaf, Zaf,'o')
xlabel('Mean Wind Speed [m/s]'); ylabel('Height [m]');
legend('Arya','Raw','Location','southeast');
xlim([0,10]);
set(gca,'FontSize',14)
% title('AFCRL - Stable')

% Plot just inlet to wind turbine
figure
plot(M_saf, z, Uaf, Zaf,'o')
xlabel('Mean Wind Speed [m/s]'); ylabel('Height [m]');
legend('Arya','Raw','Location','southeast');
% xlim([0,10]);
set(gca,'FontSize',14)
% title('AFCRL - Stable')
ylim([hub-rad, hub+rad])



m3=(Uw(4)-Uw(3))/log(Zw(4)/Zw(3));
b = exp(log(Zw(3)) - Uw(3)/m3);

figure
semilogy(M_sw, z, Uw, Zw,'o')
xlabel('Mean Wind Speed [m/s]'); ylabel('Height [m]');
legend('Arya','Raw','Location','southeast');
xlim([0,15]);
set(gca,'FontSize',14)
% title('Wangara - Stable')

