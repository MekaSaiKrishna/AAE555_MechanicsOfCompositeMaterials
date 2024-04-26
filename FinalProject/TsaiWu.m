function f_allow = TsaiWu(strength,stress)
% Solves for N22 (tensile)

% Format: strength = [X; Xdash; Y; Ydash; S; R] % in Pa
X     = strength(1);
Xdash = strength(2);
Y     = strength(3);
Ydash = strength(4);
S     = strength(5);
R     = strength(6);

% Plane Stress Assumption (sigma33=sigma23=sigma13==0)
B = stress(1)*((1/X)-(1/Xdash)) + stress(2)*((1/Y)-(1/Ydash));
A = ((stress(1))^2)/(X*Xdash) + ((stress(2))^2)/(Y*Ydash) - ...
(stress(1)*stress(2))/(X*Xdash) + (stress(3)/S)^2;
C = -1;

syms t real

if A~=0 && B~=0
    eqn = A*t^2 + B*t + C==0;
    sol = solve(eqn,t);
else
    sol = [-Inf,Inf];
end

% Store the positive solution because N22 is in tension
f_allow = max(sol);


end