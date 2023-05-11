% Data from IGRA
% https://www1.ncdc.noaa.gov/pub/data/igra/
% USM00074494  41.6569  -69.9589   15.1 MA CHATHAM   1970 2021
%
% Model from Liu and Stevens (2022)
% https://www.pnas.org/doi/10.1073/pnas.2119369119

clear
close all

%%%%%%%% plot settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(groot,'defaultfigureposition',[400 250 750 500]);
set(0,'defaultaxesfontname','times');
cOrder = [0.00 0.00 1.00
          1.00 0.65 0.00
          0.00 0.65 0.00];
mOrder = {'o','s','^'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsamples = 3; % max 3

[raw,data] = getData('USM00074494',Nsamples);

k = 0.4; % von Karman cosntant
latitude = 41.6569;
f = pi/(21600)*sind(latitude);
z0 = 0.02;
for i=1:Nsamples
    Nt = numel(data(i).GZ);
    z = data(i).GZ(data(i).GZ<=data(i).h);
    Nz = length(z);

    xi = z/data(i).h; % Greek letter for dimensionless height
    epr = (data(i).h - data(i).zi)/2/data(i).h;
    Umag = raw(i).WSPD(1:Nz);

    Tpot = data(i).TEMP(1:Nz).*(100000./data(i).PRESS(1:Nz)).^0.286;
    for j=3:Nz-1
        if (Tpot(j+1)-Tpot(j))>=0
            break
        end
    end
    p = polyfit(z(1:j),Tpot(1:j),1);
    dTPdz = p(1);
    Tpot0 = p(2);
    N = imag(sqrt(9.80665/Tpot0*dTPdz));
    Zi = N/f;

    us = bestFricVel(f,z0,Zi,k,epr,xi,Umag);

    Ro = us/f/z0; % Rossby number

    model = getModel(Ro,Zi,k,epr,z0,xi);
    lincut = 0.7;
    WDIR0 = orientAxes(model.Und,model.Vnd,Umag/us,raw(i).WDIR(1:Nz),xi,lincut);
    data(i).U = raw(i).WSPD(1:Nt).*cosd(raw(i).WDIR(1:Nt)-WDIR0);
    data(i).V = raw(i).WSPD(1:Nt).*sind(raw(i).WDIR(1:Nt)-WDIR0);
    data(i).U = sign(mean(data(i).U(1:Nz)))*data(i).U;
    p = polyfit(z(xi<lincut),data(i).V(xi<lincut),1);
    data(i).V = -sign(p(1))*data(i).V;
    Ugnd = mean(data(i).U(data(i).GZ>=data(i).h))/us;
    Vgnd = mean(data(i).V(data(i).GZ>=data(i).h))/us;

    model = getModel(Ro,Zi,k,epr,z0);

    %
    figure(3)

    subplot(1,2,1)
    hold on
    scatter(Umag/us,xi,[],cOrder(i,:),mOrder{i},'filled')
    plot(sqrt(model.Und.^2+model.Vnd.^2),model.xi,'Color',cOrder(i,:))
    pbaspect([3,4,1])
    xlabel('$U_{mag}/u_*$','FontSize',14,'Interpreter','latex')
    ylabel('$z/h$','FontSize',14,'Interpreter','latex')

    subplot(1,2,2)
    hold on
    scatter(data(i).TEMP(1:Nz),xi,[],cOrder(i,:),mOrder{i},'filled')
    plot(xlim,[0,0]+data(i).zi/data(i).h,'Color',cOrder(i,:))
    pbaspect([3,4,1])
    xlabel('$T$ (C)','FontSize',14,'Interpreter','latex')
    ylabel('$z/h$','FontSize',14,'Interpreter','latex')
    %

    %
    figure(4)
    subplot(1,2,1)
    hold on
    scatter(data(i).U(1:Nz)/us,xi,[],cOrder(i,:),mOrder{i},'filled')
    plot(model.Und,model.xi,'Color',cOrder(i,:))
    pbaspect([3,4,1])
    xlabel('$U/u_*$','FontSize',14,'Interpreter','latex')
    ylabel('$z/h$','FontSize',14,'Interpreter','latex')

    subplot(1,2,2)
    hold on
    scatter(data(i).V(1:Nz)/us,xi,[],cOrder(i,:),mOrder{i},'filled')
    plot(model.Vnd,model.xi,'Color',cOrder(i,:))
    pbaspect([3,4,1])
    xlabel('$V/u_*$','FontSize',14,'Interpreter','latex')
    ylabel('$z/h$','FontSize',14,'Interpreter','latex')
    %

    %%%%%%%% Modification of the BL height to see if it fits better %%%%%%%%
    h = data(i).zi;
    z = data(i).GZ(data(i).GZ<=h);
    Nz = length(z);
    xi = z/h;
    Umag = raw(i).WSPD(1:Nz);
    Zi = 51 + (Zi - 1245.3)/(1339.9 - 1245.3)*(89 - 51);
    us = bestFricVel(f,z0,Zi,k,epr,xi,Umag);
    Ro = us/f/z0;
    model = getModel(Ro,Zi,k,epr,z0,xi);
    lincut = 0.9;
    WDIR0 = orientAxes(model.Und,model.Vnd,Umag/us,raw(i).WDIR(1:Nz),xi,lincut);
    data(i).U = raw(i).WSPD(1:Nt).*cosd(raw(i).WDIR(1:Nt)-WDIR0);
    data(i).V = raw(i).WSPD(1:Nt).*sind(raw(i).WDIR(1:Nt)-WDIR0);
    data(i).U = sign(mean(data(i).U(1:Nz)))*data(i).U;
    p = polyfit(z(xi<lincut),data(i).V(xi<lincut),1);
    data(i).V = -sign(p(1))*data(i).V;
    Ugnd = mean(data(i).U(data(i).GZ>=data(i).h));
    Vgnd = mean(data(i).V(data(i).GZ>=data(i).h));
    model = getModel(Ro,Zi,k,epr,z0);

    %
    figure(5)

    subplot(1,2,1)
    hold on
    scatter(Umag/us,xi,[],cOrder(i,:),mOrder{i},'filled')
    plot(sqrt(model.Und.^2+model.Vnd.^2),model.xi,'Color',cOrder(i,:))
    pbaspect([3,4,1])
    xlabel('$U_{mag}/u_*$','FontSize',14,'Interpreter','latex')
    ylabel('$z/h$','FontSize',14,'Interpreter','latex')

    subplot(1,2,2)
    hold on
    scatter(data(i).TEMP(1:Nz),xi,[],cOrder(i,:),mOrder{i},'filled')
    pbaspect([3,4,1])
    xlabel('$T$ (C)','FontSize',14,'Interpreter','latex')
    ylabel('$z/h$','FontSize',14,'Interpreter','latex')
    %

    %
    figure(6)

    subplot(1,2,1)
    hold on
    scatter(data(i).U(1:Nz)/us,xi,[],cOrder(i,:),mOrder{i},'filled')
    plot(model.Und,model.xi,'Color',cOrder(i,:))
    pbaspect([3,4,1])
    xlabel('$U/u_*$','FontSize',14,'Interpreter','latex')
    ylabel('$z/h$','FontSize',14,'Interpreter','latex')

    subplot(1,2,2)
    hold on
    scatter(data(i).V(1:Nz)/us,xi,[],cOrder(i,:),mOrder{i},'filled')
    plot(model.Vnd,model.xi,'Color',cOrder(i,:))
    pbaspect([3,4,1])
    xlabel('$V/u_*$','FontSize',14,'Interpreter','latex')
    ylabel('$z/h$','FontSize',14,'Interpreter','latex')
    %
end

figure(1)
exportgraphics(gcf,[pwd '/Figures/TropopauseMethod.eps'])

figure(3)
exportgraphics(gcf,[pwd '/Figures/Umag_BL.eps'])

figure(4)
exportgraphics(gcf,[pwd '/Figures/U_V_BL.eps'])

figure(5)
exportgraphics(gcf,[pwd '/Figures/Umag_BLmod.eps'])

figure(6)
exportgraphics(gcf,[pwd '/Figures/U_V_BLmod.eps'])
