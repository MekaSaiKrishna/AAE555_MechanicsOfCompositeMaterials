% Question 5.6

% Layup: [+30/-30]s => [+30/-30/-30/+30]

% Lamina Input Properties:
E1   = 140.0e9; %in Pa = N/m^2
E2   = 10.0e9;  %in Pa = N/m^2
% E3   = E2; 
G12  = 7.0e9;   %in Pa = N/m^2
% G13  = G12;
nu12 = 0.3;
% nu13 = nu12;
% nu23 = 0.4;

% G23 = E2*0.5/(1+nu23);%in Pa

alf11 = -1e-6; % per degree centigrade
alf22 = 26e-6; % per degree centigrade

t = 0.127*1e-3; %in m (thickness of ply)

del = -100; % Temperature drop (in degrees centigrade)
%==================================================================
% Case-1: Fixed on all sides

% Find the reaction forces and laminar stresses.



%==================================================================
% Case-2: Free on all sides

% Find the plate strains and laminar stresses.



















