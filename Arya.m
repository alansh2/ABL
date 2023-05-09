% Arya Stable BL approximation
function [ustar_s, L_s, M_s]  = Arya(z1, u1, u2, Th1, Th2, k, z0, z)
    Del_u = u2-u1;
    Del_Theta = Th2-Th1;
    Theta = Th1;
    g = 9.81;
    
    ustar_s = (k*u1 - (4.7*k*g*Del_Theta*z1/(Theta*Del_u)))/log(z1/z0);
    L_s = ustar_s * Theta * Del_u/(k*g*Del_Theta);
    
    
    
    Psi_s = 4.7.*z ./L_s;
    
    M_s = ustar_s.*(1/k) .* (log(z./z0) + Psi_s);
end