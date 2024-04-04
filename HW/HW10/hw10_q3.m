% Homework-10: Question 5.9

% Layup: Antisymmetric: [60/30/-30/-60]

% Input Properties:
E1   = 140e9; %in Pa = N/m^2
E2   = 10e9;  %in Pa = N/m^2
G12  = 7e9;   %in Pa = N/m^2
nu12 = 0.3;

nu21 = E2*(nu12/E1);

t = 0.127*1e-3; %in m (thickness of ply)

alf11 = 10e-7; % Units: (/F)
alf22 = 20e-6; % Units: (/F)

%-------------------------------------------------------------------
% Step1: Find [Q] for material in the material CSYS
Q11 = E1/(1-nu12*nu21);
Q22 = E2/(1-nu12*nu21);
Q12 = nu12*E2/(1-nu12*nu21);
Q66 = G12;

Q = [Q11, Q12, 0;
    Q12, Q22, 0;
    0, 0, Q66]; %in Pa

%-------------------------------------------------------------------
% Step-2: Find Qbar for different angles
[Qd_p30] = findQbar(Q,30);  %in Pa = N/m^2
[Qd_p60] = findQbar(Q,60);  %in Pa = N/m^2
[Qd_m30] = findQbar(Q,-30); %in Pa = N/m^2
[Qd_m60] = findQbar(Q,-60); %in Pa = N/m^2

%-------------------------------------------------------------------
% Step-3: Find [A], [B], [D] matrices
syms x real
z0 = -2*t; z1 = -t; z2 = 0; 
z3 = t; z4 = 2*t; %in m 


Q1 = Qd_p60;  Q2 = Qd_p30;   %in Pa
Q3 = Qd_m30;  Q4 = Qd_m60;   %in Pa


[A] = vpa(int(Q1,x,z0,z1),4) + vpa(int(Q2,x,z1,z2),4) + vpa(int(Q3,x,z2,z3),4)...
    + vpa(int(Q4,x,z3,z4),4);  %in = (N/m^2)*m = (N/m)

A = vpa(A,4) %in = (N/m^2)*m = (N/m)

[B] = vpa(int(Q1*x,x,z0,z1) + int(Q2*x,x,z1,z2) + int(Q3*x,x,z2,z3) ...
    + int(Q4*x,x,z3,z4),4); %in Pa*m^2 = (N/m^2)*(m)^2  = N

B = vpa(B,4) %in Pa*m^2 = (N/m^2)*(m)^2  = N

[D] = vpa(int(Q1*x^2,x,z0,z1) + int(Q2*x^2,x,z1,z2) + int(Q3*x^2,x,z2,z3) ...
    + int(Q4*x^2,x,z3,z4),4); %in Pa*m^3 = (N/m^2)*(m)^3 = (N*m)

D = vpa(D,4) %in Pa*m^3 = (N/m^2)*(m)^3 = (N*m)


%-------------------------------------------------------------------
% Step-4: Find [NTbar], [MTbar] vectors

% Step-4.0: Thermal Coefficients in-plane
alf0 = [alf11,alf22, 0]';

% Step-4.1: Find alpha-dash for different angles
R_p60 = R_eps_e(+60); R_p30 = R_eps_e(+30);
R_m60 = R_eps_e(-60); R_m30 = R_eps_e(-30);

alf_p60 = R_p60*alf0; alf_m60 = R_m60*alf0;
alf_p30 = R_p30*alf0; alf_m30 = R_m30*alf0;

% Assign alphaDash for each lamina
alf1 = alf_p60; alf2 = alf_p30;
alf3 = alf_m30; alf4 = alf_m60;

NTbar = (vpa(int(Q1*alf1,x,z0,z1),4) + vpa(int(Q2*alf2,x,z1,z2),4) + ...
    vpa(int(Q3*alf3,x,z2,z3),4) + vpa(int(Q4*alf4,x,z3,z4),4)); 

NTbar=vpa(NTbar,6) %Units: (N/(m*F))

MTbar = (vpa(int(Q1*alf1*x,x,z0,z1),4) + vpa(int(Q2*alf2*x,x,z1,z2),4) + ...
    vpa(int(Q3*alf3*x,x,z2,z3),4) + vpa(int(Q4*alf4*x,x,z3,z4),4));

MTbar=vpa(MTbar,6) %Units: (N/F)

%-------------------------------------------------------------------
% Step-5: Find [NTbar_dash], [MTbar_dash] vectors

NTbar_dash = (vpa(int(Q1*alf1*x,x,z0,z1),4) + vpa(int(Q2*alf2*x,x,z1,z2),4) + ...
    vpa(int(Q3*alf3*x,x,z2,z3),4) + vpa(int(Q4*alf4*x,x,z3,z4),4));

NTbar_dash=vpa(NTbar_dash,6) %Units: (N/F)

MTbar_dash = (vpa(int(Q1*alf1*x^2,x,z0,z1),4) + vpa(int(Q2*alf2*x^2,x,z1,z2),4) + ...
    vpa(int(Q3*alf3*x^2,x,z2,z3),4) + vpa(int(Q4*alf4*x^2,x,z3,z4),4));

MTbar_dash=vpa(MTbar_dash,6) %Units: (N*m/F)



















