% Homework-10: Question 5.8

% Layup: [+45/-45/0]s => [+45/-45/0/0/-45/+45]

% Input Properties:
E1   = 140e9; %in Pa = N/m^2
E2   = 10e9;  %in Pa = N/m^2
G12  = 7e9;   %in Pa = N/m^2
nu12 = 0.3;

nu21 = E2*(nu12/E1);

t = 0.127*1e-3; %in m (thickness of ply)

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
[Qd_0]   = findQbar(Q,0);   %in Pa = N/m^2
[Qd_p45] = findQbar(Q,45);  %in Pa = N/m^2
[Qd_m45] = findQbar(Q,-45); %in Pa = N/m^2

%-------------------------------------------------------------------
% Step-3: Find [D] matrix
syms x real
z0 = -3*t; z1 = -2*t; z2 = -t; z3 = 0; %in m
z4 = t;    z5 = 2*t;  z6 = 3*t; %in m 

%Layup: [0/60/-60]s => [0/60/-60/-60/60/0]

Q1 = Qd_p45;  Q2 = Qd_m45;  Q3 = Qd_0;     %in Pa
Q4 = Qd_0;    Q5 = Qd_m45;  Q6 = Qd_p45;   %in Pa


[D] = vpa(int(Q1*x^2,x,z0,z1),4) + vpa(int(Q2*x^2,x,z1,z2),4) + vpa(int(Q3*x^2,x,z2,z3),4)...
    + vpa(int(Q4*x^2,x,z3,z4),4) + vpa(int(Q5*x^2,x,z4,z5),4) + ...
      vpa(int(Q6*x^2,x,z5,z6),4); %in Pa*m^3 = (N/m^2)*(m)^3 = (N*m)

D = vpa(D,4); %in Pa*m^3 = (N/m^2)*(m)^3 = (N*m)

%--------------------------------------------------------------------------
% Step-4: Find (flexural) Effective Plane Stress Reduced Compliance Matrix 
%         [Se_starf] 

h = (6*t);

Se_starf = (h^3/12)*inv(D); %in m^2/N

%--------------------------------------------------------------------------
% Step-5: Effective In-plane Properties (flexural)
[E1_flex,E2_flex,nu12_flex,G12_flex,eta121_flex,eta122_flex] = findEffInPlane(Se_starf);


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fprintf("      Effective In-Plane Properties     \n")
fprintf(" ====================================== \n")
fprintf(" Prop     |      CLT   (flexural)       \n")
fprintf(" -------------------------------------- \n")
fprintf(" E1       | %0.2f  GPa                  \n", vpa(E1_flex*1e-9,2))
fprintf(" E2       | %0.2f  GPa                  \n", vpa(E2_flex*1e-9,2))
fprintf(" nu12     | %0.2f  GPa                  \n", vpa(nu12_flex,2))
fprintf(" G12      | %0.2f  GPa                  \n", vpa(G12_flex*1e-9,2))
fprintf(" eta121   | %0.2f                       \n", vpa(eta121_flex,2))
fprintf(" eta122   | %0.2f                       \n", vpa(eta122_flex,2))





















