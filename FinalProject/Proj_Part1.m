% Part-1: HROM for Lamina Properties

% Fiber (Carbon) Properties
global E1 E2 nu12 nu23 G12
E1   = 276.0*1e9; % in Pa
E2   = 19.5*1e9;  % in Pa
nu12 = 0.28;
nu23 = 0.70;
G12  = 70.0*1e9; % in Pa


% Matrix (Epoxy) Properties
global E nu 

E  = 4.76*1e9; % in Pa
nu = 0.37;

% Volume fraction of fiber
Vf = [0.5,0.6,0.7];

% Calculate Effective Engineering constants
% using HROM for different Volume Fractions

% Store the HROM engineering constants for all volume fractions
CONSTS = zeros(5,3);

for i=1:3
    PROPS = HROM(Vf(i));
    fprintf("Fiber Volume Fraction: %0.2f \n",Vf(i))
    fprintf("------------------------------------- \n")
    fprintf("E1H*   : %0.2f GPa \n", 1e-9*PROPS(1))
    fprintf("nu12H* : %0.2f \n", PROPS(2))
    fprintf("E2H*   : %0.2f GPa \n", 1e-9*PROPS(3))
    fprintf("nu23H* : %0.2f \n", PROPS(4))
    fprintf("G12H*  : %0.2f GPa \n", 1e-9*PROPS(5))
    fprintf("===================================== \n")
    fprintf("\n")

    CONSTS(:,i) = PROPS;
end

save("consts.mat","CONSTS")



















