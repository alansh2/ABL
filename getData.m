function [raw,data] = getData(station,Nsamples)

    raw = struct('GPH',cell(Nsamples,1),...
        'PRESS',cell(Nsamples,1),...
        'TEMP',cell(Nsamples,1),...
        'WDIR',cell(Nsamples,1),...
        'WSPD',cell(Nsamples,1));
    data = struct('GZ',cell(Nsamples,1),...
        'PRESS',cell(Nsamples,1),...
        'TEMP',cell(Nsamples,1),...
        'U',cell(Nsamples,1),...
        'V',cell(Nsamples,1),...
        'h',cell(Nsamples,1),...
        'zi',cell(Nsamples,1));

    getRaw;

    for I=1:Nsamples
        gzh = gp2gz(raw(I).GPH);
        J = findTropopause(gzh,raw(I).TEMP);
        % Add to data the variables that do not need further processing
        data(I).GZ = gzh(J);
        data(I).TEMP = raw(I).TEMP(J);
        data(I).PRESS = raw(I).PRESS(J);
        % Compute wind variables
        data(I).U = raw(I).WSPD(J).*cosd(raw(I).WDIR(J));
        data(I).V = raw(I).WSPD(J).*sind(raw(I).WDIR(J));

        [data(I).h,data(I).zi] = getBndLyrHgts(data(I).GZ,data(I).TEMP);
    end

    fig = figure(1);
    han = axes(fig,'visible','off');
    % Matlab only
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    %
    ylabel(han,'Geometric Height (m)','FontSize',14);
    xlabel(han,'Temperature (C)','FontSize',14);

    %% SUBROUTINES
    function getRaw
        fid = fopen(['Data\' station '-data-beg2021.txt']);
        K = 1;
        while ~feof(fid)
            tmp = fgetl(fid);
            if (tmp(1)=='#')
                if (size(raw(K).GPH,1)<30)
                    raw(K).PRESS = [];
                    raw(K).GPH   = [];
                    raw(K).TEMP  = [];
                    raw(K).WDIR  = [];
                    raw(K).WSPD  = [];
                else
                    K = K + 1;
                end
                continue
            end

            if (K>Nsamples)
                break
            end

            raw(K).PRESS = [raw(K).PRESS;str2double(tmp(10:15))];
            raw(K).GPH   = [raw(K).GPH;str2double(tmp(17:21))];
            raw(K).TEMP  = [raw(K).TEMP;str2double(tmp(23:27))/10];
            raw(K).WDIR  = [raw(K).WDIR;str2double(tmp(41:45))];
            raw(K).WSPD  = [raw(K).WSPD;str2double(tmp(47:51))/10];
        end
        fclose(fid);
    end

    function gzhgt = gp2gz(gphgt)
        % Convert geopotential height to geometric height
        G  = 6.674e-11; % gravitational constant
        M  = 5.9722e24; % mass of Earth
        R  = 6.3781e6;  % radius of Earth
        g0 = 9.80665;   % reference gravitational acceleration

        gpot = gphgt*g0; % geopotential
        gzhgt = 1./(1/R-gpot/(G*M)) - R; % geometric height
    end

    function ndx = findTropopause(gz,T)
        J = find((gz>7000)&(gz<20000));
        N = J(end) - J(1) + 1;

        % Fit a sigmoid function to the temperature profile
        lb = [20 -9e-3 7000 -40];
        ub = [75 -6e-3 20000 0];
        x0 = [35 -6.5e-3 9000 -30];
        x0 = scaling(x0,lb,ub,1);
        xopt = fmincon(@fun,x0,[],[],[],[],zeros(1,4),ones(1,4));
%        xopt = fminunc(@fun,x0); % unconstrained optimization possible
        xopt = scaling(xopt,lb,ub,2);

        L = xopt(2); % temperature lapse rate
        Tmin = xopt(4) - xopt(1);
        hTrop = xopt(3) - (xopt(4)-Tmin)/L; % height where Tropopause begins

        ndx = gz<hTrop;

        figure(1)
        subplot(1,Nsamples,I);
        hold on
        scatter(T(gz<20000),gz(gz<20000),[],'blue')
        plot(xopt(1)*tanh(xopt(2)/xopt(1)*((7000:100:20000)-xopt(3)))+xopt(4),(7000:100:20000),'k--')
        plot([Tmin Tmin -20],[20000 hTrop xopt(3)-(xopt(4)+20)/L],'k-')
        xlim([-70,10])

        function y = fun(xin)
            x = scaling(xin,lb,ub,2);
            y = sqrt(1/N*sum((x(1)*tanh(x(2)/x(1)*(gz(J)-x(3)))+x(4)-T(J)).^2));
        end

    end

    function x_scaled = scaling(x,LB,UB,s_type)
        if s_type==1
            x_scaled = (x - LB)./(UB - LB);
        elseif s_type==2
            x_scaled = LB + x.*(UB - LB);
        else
            msg  = 'Invalid scaling operation specified. Input 1 to scale or 2 to unscale.';
            error(msg)
        end
    end

    function [h,zi] = getBndLyrHgts(gz,T)
        idx = find(gz<2700);
        N = idx(end) - idx(1) + 1;
        gz11 = 2*scaling(gz(idx),0,2700,1) - 1; % [-1,1]
        ugz = linspace(-1,1,N)';
        pdeg = 7;
        p = polyfit(ugz,T(idx),pdeg);
        d = p(1:pdeg).*(pdeg:-1:1);
        dd = d(1:pdeg-1).*(pdeg-1:-1:1);
        r = roots(d);
        % Remove imaginary and OOB roots and sort in increasing height
        r = r(imag(r)==0);
        r = r((r>=-1)&(r<=1));
        r = sort(r);
        % Find roots with upward concavity
        ndx = find(polyval(dd,r)>0);
        I1 = ndx(1); % first guess at zi (unused)
        I2 = I1 + 1; % for h
        A = sort([ugz;r(I1)]);
        B = sort([ugz;r(I2)]);
%        zi = scaling((gz11(find(A==r(I1)))+1)/2,0,2700,2);
        h = scaling((gz11(B==r(I2))+1)/2,0,2700,2);

        %%%%%%%% Refine the estimate of zi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx = find(gz<=h);
        N = idx(end) - idx(1) + 1;
        gz11 = 2*scaling(gz(idx),0,gz(idx(end)),1) - 1; % [-1,1]
        pdeg = 6;
        p = polyfit(gz11,T(idx),pdeg);
        d = p(1:pdeg).*(pdeg:-1:1);
        dd = d(1:pdeg-1).*(pdeg-1:-1:1);
        r = roots(d);
        % Remove imaginary and OOB roots and sort in increasing height
        r = r(imag(r)==0);
        r = r((r>=-1)&(r<=1));
        r = sort(r);
        % Find roots with upward concavity
        ndx = find(polyval(dd,r)>0);
        I1 = ndx(1); % refined guess of zi
        zi = scaling((r(I1)+1)/2,0,gz(idx(end)),2);
%        figure(2)
%        hold on
%        scatter(T(idx),gz11)
%        plot(polyval(p,linspace(-1,1,201)),linspace(-1,1,201))
    end

end
