% Compute effective engineering constants using
% HROM 

function PROPS = HROM(Vf)

% Input: 
% Vf: Fiber volume fraction

% Output:
% PROPS(1) = E1H*;
% PROPS(2) = nu12H*;
% PROPS(3) = E2H*;
% PROPS(4) = nu23H*;
% PROPS(5) = G12H*;

global E1 E2 nu12 nu23 G12
global E nu 

Gm = 0.5*E/(1+nu);

Vm = 1-Vf;

E1H   = Vf*E1 + Vm*E;
nu12H = Vf*nu12 + Vm*nu;

%-----------------------------------
term1 = Vf/E2; term2 = Vm/E;
term3 = (Vf*Vm*(E*nu12 - E1*nu)^2)/(E1*E*(E1*Vf + E*Vm));

E2H = 1/(term1+term2-term3);
%-----------------------------------
term4 = Vf*((nu23/E2)+((nu12^2)/E1));
term5 = Vm*nu*(1+nu)/E;
term6 = nu12H^2/E1H;

nu23H = E2H*(term4+term5-term6);

%-----------------------------------
G12H = 1/((Vf/G12)+(Vm/Gm));

PROPS = [E1H, nu12H, E2H, nu23H, G12H];

end
