% Homework-10: Question 5.10

% Find possible stack sequences to satisfy Quasi-isotropic requirement
% Conditions:
%   (1) 6 layer laminate
%   (2) Symmetric (Because Quasi Isotropic)

%-------------------------------------------------------------------
% Possible Combinations (3!):
% Combination-1: [0,60,-60]s = [0/60/-60/-60/60/0]
% Combination-2: [0,-60,60]s = [0/-60/60/60/-60/0]
% Combination-3: [60,0,-60]s = [60/0/-60/-60/0/60]
% Combination-4: [-60,0,60]s = [-60/0/60/60/0/-60]
% Combination-5: [60,-60,0]s = [60/-60/0/0/-60/60]
% Combination-6: [-60,60,0]s = [-60/60/0/0/60/-60]


%-------------------------------------------------------------------

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
[Qd_0]   = findQbar(Q,0);  %in Pa = N/m^2
[Qd_p60] = findQbar(Q,60);  %in Pa = N/m^2
[Qd_m60] = findQbar(Q,-60); %in Pa = N/m^2


% Possible Combinations (3!):

%====================================================
% Combination-1: [0,60,-60]s = [0/60/-60/-60/60/0]

% Step-3: Find [A], [D] matrices
Q1 = Qd_0; Q2 = Qd_p60; Q3 = Qd_m60;
[A1,D1] = findAD(Q1,Q2,Q3,t);

%--------------------------------------------------------------------------
% Step-4: Find Effective Plane Stress Reduced Compliance Matrix [Se_star] 

h = (6*t);
Se_star = vpa(h*inv(A1)); %in m^2/N
Se_starf = (h^3/12)*inv(D1); %in m^2/N

%--------------------------------------------------------------------------
% Step-5: Effective In-plane Properties
[E1_star,E2_star,nu12_star,G12_star,eta121_star,eta122_star] = findEffInPlane(Se_star);
[E1_flex,E2_flex,nu12_flex,G12_flex,eta121_flex,eta122_flex] = findEffInPlane(Se_starf);

fprintf("      Combination-1: [0,60,-60]s       \n")
fprintf(" ====================================== \n")
fprintf(" Prop     |    INPLANE  |      FLEX     \n")
fprintf(" -------------------------------------- \n")
fprintf(" E1       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E1_star*1e-9,2), vpa(E1_flex*1e-9,2))
fprintf(" E2       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E2_star*1e-9,2), vpa(E2_flex*1e-9,2))
fprintf(" nu12     | %0.2f        | %0.2f        \n", vpa(nu12_star,2), vpa(nu12_flex,2))
fprintf(" G12      | %0.2f  GPa  | %0.2f GPa     \n", vpa(G12_star*1e-9,2), vpa(G12_flex*1e-9,2))
fprintf(" eta121   | %0.2f        | %0.2f        \n", vpa(eta121_star,2), vpa(eta121_flex,2))
fprintf(" eta122   | %0.2f        | %0.2f        \n", vpa(eta122_star,2), vpa(eta122_flex,2))


%====================================================
% Combination-2: [0,-60,60]s = [0/-60/60/60/-60/0]

% Step-3: Find [A], [D] matrices
Q1 = Qd_0; Q2 = Qd_m60; Q3 = Qd_p60;
[A2,D2] = findAD(Q1,Q2,Q3,t);

%--------------------------------------------------------------------------
% Step-4: Find Effective Plane Stress Reduced Compliance Matrix [Se_star] 

h = (6*t);
Se_star = vpa(h*inv(A2)); %in m^2/N
Se_starf = (h^3/12)*inv(D2); %in m^2/N

%--------------------------------------------------------------------------
% Step-5: Effective In-plane Properties
[E1_star,E2_star,nu12_star,G12_star,eta121_star,eta122_star] = findEffInPlane(Se_star);
[E1_flex,E2_flex,nu12_flex,G12_flex,eta121_flex,eta122_flex] = findEffInPlane(Se_starf);

fprintf(" \n")
fprintf(" \n")
fprintf("      Combination-2: [0,-60,60]s       \n")
fprintf(" ====================================== \n")
fprintf(" Prop     |    INPLANE  |      FLEX     \n")
fprintf(" -------------------------------------- \n")
fprintf(" E1       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E1_star*1e-9,2), vpa(E1_flex*1e-9,2))
fprintf(" E2       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E2_star*1e-9,2), vpa(E2_flex*1e-9,2))
fprintf(" nu12     | %0.2f        | %0.2f        \n", vpa(nu12_star,2), vpa(nu12_flex,2))
fprintf(" G12      | %0.2f  GPa  | %0.2f GPa     \n", vpa(G12_star*1e-9,2), vpa(G12_flex*1e-9,2))
fprintf(" eta121   | %0.2f        | %0.2f        \n", vpa(eta121_star,2), vpa(eta121_flex,2))
fprintf(" eta122   | %0.2f        | %0.2f        \n", vpa(eta122_star,2), vpa(eta122_flex,2))


%====================================================
% Combination-3: [60,0,-60]s = [60/0/-60/-60/0/60]

% Step-3: Find [A], [D] matrices
Q1 = Qd_p60; Q2 = Qd_0; Q3 = Qd_m60;
[A3,D3] = findAD(Q1,Q2,Q3,t);

%--------------------------------------------------------------------------
% Step-4: Find Effective Plane Stress Reduced Compliance Matrix [Se_star] 

h = (6*t);
Se_star = vpa(h*inv(A3)); %in m^2/N
Se_starf = (h^3/12)*inv(D3); %in m^2/N

%--------------------------------------------------------------------------
% Step-5: Effective In-plane Properties
[E1_star,E2_star,nu12_star,G12_star,eta121_star,eta122_star] = findEffInPlane(Se_star);
[E1_flex,E2_flex,nu12_flex,G12_flex,eta121_flex,eta122_flex] = findEffInPlane(Se_starf);

fprintf(" \n")
fprintf(" \n")
fprintf("      Combination-3: [60,0,-60]s       \n")
fprintf(" ====================================== \n")
fprintf(" Prop     |    INPLANE  |      FLEX     \n")
fprintf(" -------------------------------------- \n")
fprintf(" E1       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E1_star*1e-9,2), vpa(E1_flex*1e-9,2))
fprintf(" E2       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E2_star*1e-9,2), vpa(E2_flex*1e-9,2))
fprintf(" nu12     | %0.2f        | %0.2f        \n", vpa(nu12_star,2), vpa(nu12_flex,2))
fprintf(" G12      | %0.2f  GPa  | %0.2f GPa     \n", vpa(G12_star*1e-9,2), vpa(G12_flex*1e-9,2))
fprintf(" eta121   | %0.2f        | %0.2f        \n", vpa(eta121_star,2), vpa(eta121_flex,2))
fprintf(" eta122   | %0.2f        | %0.2f        \n", vpa(eta122_star,2), vpa(eta122_flex,2))


%====================================================
% Combination-4: [-60,0,60]s = [-60/0/60/60/0/-60]

% Step-3: Find [A], [D] matrices
Q1 = Qd_m60; Q2 = Qd_0; Q3 = Qd_p60;
[A4,D4] = findAD(Q1,Q2,Q3,t);

%--------------------------------------------------------------------------
% Step-4: Find Effective Plane Stress Reduced Compliance Matrix [Se_star] 

h = (6*t);
Se_star = vpa(h*inv(A4)); %in m^2/N
Se_starf = (h^3/12)*inv(D4); %in m^2/N

%--------------------------------------------------------------------------
% Step-5: Effective In-plane Properties
[E1_star,E2_star,nu12_star,G12_star,eta121_star,eta122_star] = findEffInPlane(Se_star);
[E1_flex,E2_flex,nu12_flex,G12_flex,eta121_flex,eta122_flex] = findEffInPlane(Se_starf);

fprintf(" \n")
fprintf(" \n")
fprintf("      Combination-4: [-60,0,60]s       \n")
fprintf(" ====================================== \n")
fprintf(" Prop     |    INPLANE  |      FLEX     \n")
fprintf(" -------------------------------------- \n")
fprintf(" E1       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E1_star*1e-9,2), vpa(E1_flex*1e-9,2))
fprintf(" E2       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E2_star*1e-9,2), vpa(E2_flex*1e-9,2))
fprintf(" nu12     | %0.2f        | %0.2f        \n", vpa(nu12_star,2), vpa(nu12_flex,2))
fprintf(" G12      | %0.2f  GPa  | %0.2f GPa     \n", vpa(G12_star*1e-9,2), vpa(G12_flex*1e-9,2))
fprintf(" eta121   | %0.2f        | %0.2f        \n", vpa(eta121_star,2), vpa(eta121_flex,2))
fprintf(" eta122   | %0.2f        | %0.2f        \n", vpa(eta122_star,2), vpa(eta122_flex,2))

%====================================================
% Combination-5: [60,-60,0]s = [60/-60/0/0/-60/60]

% Step-3: Find [A], [D] matrices
Q1 = Qd_p60; Q2 = Qd_m60; Q3 = Qd_0;
[A5,D5] = findAD(Q1,Q2,Q3,t);

%--------------------------------------------------------------------------
% Step-4: Find Effective Plane Stress Reduced Compliance Matrix [Se_star] 

h = (6*t);
Se_star = vpa(h*inv(A5)); %in m^2/N
Se_starf = (h^3/12)*inv(D5); %in m^2/N

%--------------------------------------------------------------------------
% Step-5: Effective In-plane Properties
[E1_star,E2_star,nu12_star,G12_star,eta121_star,eta122_star] = findEffInPlane(Se_star);
[E1_flex,E2_flex,nu12_flex,G12_flex,eta121_flex,eta122_flex] = findEffInPlane(Se_starf);

fprintf(" \n")
fprintf(" \n")
fprintf("      Combination-5: [60,-60,0]s       \n")
fprintf(" ====================================== \n")
fprintf(" Prop     |    INPLANE  |      FLEX     \n")
fprintf(" -------------------------------------- \n")
fprintf(" E1       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E1_star*1e-9,2), vpa(E1_flex*1e-9,2))
fprintf(" E2       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E2_star*1e-9,2), vpa(E2_flex*1e-9,2))
fprintf(" nu12     | %0.2f        | %0.2f        \n", vpa(nu12_star,2), vpa(nu12_flex,2))
fprintf(" G12      | %0.2f  GPa  | %0.2f GPa     \n", vpa(G12_star*1e-9,2), vpa(G12_flex*1e-9,2))
fprintf(" eta121   | %0.2f        | %0.2f        \n", vpa(eta121_star,2), vpa(eta121_flex,2))
fprintf(" eta122   | %0.2f        | %0.2f        \n", vpa(eta122_star,2), vpa(eta122_flex,2))


%====================================================
% Combination-6: [-60,60,0]s = [-60/60/0/0/60/-60]

% Step-3: Find [A], [D] matrices
Q1 = Qd_m60; Q2 = Qd_p60; Q3 = Qd_0;
[A6,D6] = findAD(Q1,Q2,Q3,t);

%--------------------------------------------------------------------------
% Step-4: Find Effective Plane Stress Reduced Compliance Matrix [Se_star] 

h = (6*t);
Se_star = vpa(h*inv(A6)); %in m^2/N
Se_starf = (h^3/12)*inv(D6); %in m^2/N

%--------------------------------------------------------------------------
% Step-5: Effective In-plane Properties
[E1_star,E2_star,nu12_star,G12_star,eta121_star,eta122_star] = findEffInPlane(Se_star);
[E1_flex,E2_flex,nu12_flex,G12_flex,eta121_flex,eta122_flex] = findEffInPlane(Se_starf);

fprintf(" \n")
fprintf(" \n")
fprintf("      Combination-6: [-60,60,0]s       \n")
fprintf(" ====================================== \n")
fprintf(" Prop     |    INPLANE  |      FLEX     \n")
fprintf(" -------------------------------------- \n")
fprintf(" E1       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E1_star*1e-9,2), vpa(E1_flex*1e-9,2))
fprintf(" E2       | %0.2f  GPa  | %0.2f GPa     \n", vpa(E2_star*1e-9,2), vpa(E2_flex*1e-9,2))
fprintf(" nu12     | %0.2f        | %0.2f        \n", vpa(nu12_star,2), vpa(nu12_flex,2))
fprintf(" G12      | %0.2f  GPa  | %0.2f GPa     \n", vpa(G12_star*1e-9,2), vpa(G12_flex*1e-9,2))
fprintf(" eta121   | %0.2f        | %0.2f        \n", vpa(eta121_star,2), vpa(eta121_flex,2))
fprintf(" eta122   | %0.2f        | %0.2f        \n", vpa(eta122_star,2), vpa(eta122_flex,2))




