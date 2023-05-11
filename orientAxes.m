function WDIR0 = orientAxes(Uref,Vref,Umag,Wdir,xi,lincut)
    if nargin == 5
        lincut = 0.8;
    end

    WDIR0 = fminunc(@objfun,0);

    function y = objfun(x)
        U = Umag.*cosd(Wdir - x);
        U = sign(mean(U))*U;
        V = Umag.*sind(Wdir - x);
        p = polyfit(xi(xi<lincut),V(xi<lincut),1);
        y = abs(p(2)); % cost function to enforce zero crosswind at surface
%        V = -sign(p(1))*V;
%        y = sum(abs(U - Uref)+abs(V - Vref)); % cost function that minimmizes total error
    end

end
