function us = bestFricVel(f,z0,Zi,k,epr,xi,refDat)

    us = fmincon(@objfun,0.4,[],[],[],[],0,Inf);

    function y = objfun(x)
        model = getModel(x/f/z0,Zi,k,epr,z0,xi);
        Umag = sqrt(model.Und.^2+model.Vnd.^2);
        y = sum(abs(Umag-refDat/x));
    end

end
