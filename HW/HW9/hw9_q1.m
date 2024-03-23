% Question 5.4

N11 = 4500; %in N/m

% Input Properties:
E1   = 140e9; %in Pa = N/m^2
E2   = 10e9;  %in Pa = N/m^2
G12  = 7e9;   %in Pa = N/m^2
nu12 = 0.3;

nu21 = E2*(nu12/E1);

t    = 0.127*1e-3; %in m

% What are the corresponding moments needed to make sure that laminate 
% will not bend?

% Layup: [0/90]2 => [0/90/0/90]

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
[Qd_0]  = findQbar(Q,0);  %in Pa = N/m^2
[Qd_90] = findQbar(Q,90); %in Pa = N/m^2

%-------------------------------------------------------------------
% Step-3: Find [A, B, D] matrix
syms x real
z0 = -2*t; z1 = -t; z2 = 0; z3 = t; z4 = 2*t; %in m

%Layup1: [0/90/0/90]
[A_1] = t*(2*(Qd_0 + Qd_90)) %in = (N/m^2)*m = (N/m)

Q1 = Qd_0; Q2 = Qd_90;  Q3 = Qd_0; Q4 = Qd_90;
[B_1] = vpa(int(Q1*x,x,z0,z1) + int(Q2*x,x,z1,z2) + int(Q3*x,x,z2,z3) ...
    + int(Q4*x,x,z3,z4),4) %in Pa*m^2 = (N/m^2)*(m)^2  = N

[D_1] = vpa(int(Q1*x^2,x,z0,z1) + int(Q2*x^2,x,z1,z2) + int(Q3*x^2,x,z2,z3) ...
    + int(Q4*x^2,x,z3,z4),4) %in Pa*m^3 = (N/m^2)*(m)^3 = (N*m)

fprintf("------------------------------ \n")

%-------------------------------------------------------------------
% Step-4: Find [a, b, d] matrix
[a_1,b_1,d_1] = find_abd(A_1,B_1,D_1);

abd = [a_1, b_1;
       b_1', d_1]; % SI Units

%-------------------------------------------------------------------
% Step-5: Loading condition
syms M11 M22 M12 real

%NM_vec = [N11, N22, N12, M11, M22, M12];
NM_vec = [N11, 0, 0, M11, M22, M12]';

% What are the corresponding moments needed to make sure that laminate 
% will not bend?

% Therefore, curvatures should be zero
syms E11 E22 E12 real
%ek_vec = [E11, E22, 2*E12, k11, k22, 2*k12];
ek_vec = [E11, E22, 2*E12, 0, 0, 0]';


%-------------------------------------------------------------------
% Step-6: Plate Constitutive Relation

EQN = ek_vec - abd*NM_vec == [0,0,0,0,0,0]'


sol = solve(EQN, [E11, E22, E12, M11, M22, M12]);

E11_Sol = vpa(1e6*sol.E11,4) % micro strains
E22_Sol = vpa(1e6*sol.E22,4) % micro strains
E12_Sol = vpa(1e6*sol.E12,4) % micro strains
M11_Sol = vpa(sol.M11,4) % in N
M22_Sol = vpa(sol.M22,4) % in N
M12_Sol = vpa(sol.M12,4) % in N

fprintf("E11_sol = %0.2f microstrains \n",E11_Sol)
fprintf("E22_sol = %0.2f microstrains \n",E22_Sol)
fprintf("E12_sol = %0.2f microstrains \n",E12_Sol)
fprintf("M11_sol = %0.2f N \n",M11_Sol)
fprintf("M22_sol = %0.2f N \n",M22_Sol)
fprintf("M12_sol = %0.2f N \n",M12_Sol)

