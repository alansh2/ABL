function model = getModel(Ro,Zi,k,epr,z0,xi)
% spatial discretization
if (nargin==5)
    Nz = 100;
    xi = linspace(1/Nz,1,Nz)'; % z/h
end

% constant functions for each Zi
CR = 0.5;
CN = 1.6;
beta = sqrt(CR^(-2) + CN^(-2)*Zi);
b = (1 - 0.05^(2/3))*beta;
A1 = 0.65;
A0 = 1.3;
B1 = 7;
B0 = 8;
Cm = 0.1;
m = sqrt(1 + Cm^2*Zi.^2)./beta;
A = -A1*m + log(A0 + m) + log(beta);
B = (B0 + B1*m.^2)./beta;
a = A - log(b);

% influence of Ro
h = Ro*z0./b;
zi = h*(1 - 2*epr);

% the "psi" functions
ap = (2 - a)/(1 - 2*epr);
bp = (2*k*b - B)/(1 - 2*epr);
psi = xi - (exp(xi/epr) - 1)/(exp(1/epr) - 1);
fu = -a.*xi + ap.*psi;
fv = -B.*xi + bp.*psi;

% nondimensional geostrophic velocity
Ugnd = (log(Ro) - A)/k;
Vgnd = -B/k;

Und = (log(xi.*Ro./b) + fu)/k; % xi*Ro/b = xi*h/z0 = z/z0
Vnd = fv/k;

model = struct('xi',xi,'Ugnd',Ugnd,'Vgnd',Vgnd,'Und',Und,'Vnd',Vnd);
