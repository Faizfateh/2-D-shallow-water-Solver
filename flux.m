function [F, smag] = FluxFunction(UL, UR, n)
% PURPOSE: This routine calculates the flux for the shallow-water
% equations using the Roe flux function
%
% INPUTS:
%    UL: conservative state vector in the left cell
%    UR: conservative state vector in the right cell
%     n: unit normal pointing from the left cell to the right cell
%
% OUTPUTS:
%  F   : the flux out of the left cell (into the right cell)
%  smag: the maximum propagation speed of disturbances
%

% acceleration due to gravity
g = 9.8;

% process left state
hL = UL(1);
uL = UL(2)/hL;
vL = UL(3)/hL;
unL = uL*n(1) + vL*n(2);
pL = 0.5*g*hL^2;
if (hL<=0), error 'Non-physical state!', end

% left flux
FL = UL; % for allocation
FL(1) = hL*unL;
FL(2) = hL*uL*unL + pL*n(1);
FL(3) = hL*vL*unL + pL*n(2);

% process right state
hR = UR(1);
uR = UR(2)/hR;
vR = UR(3)/hR;
unR = uR*n(1) + vR*n(2);
pR = 0.5*g*hR^2;
if (hR<=0), error 'Non-physical state!', end

% right flux
FR = UR; % for allocation
FR(1) = hR*unR;
FR(2) = hR*uR*unR + pR*n(1);
FR(3) = hR*vR*unR + pR*n(2);

% average state
h = 0.5*(hL + hR);
hu = 0.5*(hL*uL + hR*uR);
hv = 0.5*(hL*vL + hR*vR);
u = hu/h;
v = hv/h;
un = u*n(1) + v*n(2);
c = sqrt(g*h);

% difference in states
du = UR - UL;

% eigenvalues
l(1) = un;
l(2) = un-c;
l(3) = un+c;

% entropy fix
epsilon = c*.05;
for i=1:3
  if ((l(i)<epsilon) && (l(i)>-epsilon))
    l(i) = 0.5*(epsilon + l(i)*l(i)/epsilon);
  end
end

% sign of eigs
el = sign(l);

% combination of eigenvalues
s2 = 0.5*(el(2)*l(2)-el(3)*l(3));
s3 = 0.5*(el(2)*l(2)+el(3)*l(3)-2*el(1)*l(1));

% absolute value of first eigenvalue
l1 = el(1)*l(1);

% eigenvetor product generator
G1 = du(1)*un - du(2)*n(1) - du(3)*n(2);

% functions of G1, s2, s3
C1 = du(1)*s3 + G1*s2/c;
C2 = G1*s3 + s2*du(1)*c;

% flux assembly
F = FL; % for allocation
F(1)    = 0.5*(FL(1)+FR(1))-0.5*(l1*du(1) + C1   );
F(2)    = 0.5*(FL(2)+FR(2))-0.5*(l1*du(2) + C1*u - C2*n(1));
F(3)    = 0.5*(FL(3)+FR(3))-0.5*(l1*du(3) + C1*v - C2*n(2));

% max wave speed
smag = max(abs(l));
end

