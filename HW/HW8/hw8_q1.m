%% Homework 8: AAE555

% Question 5.1:

% Input Properties:
E1   = 140.0e9; %in Pa
E2   = 10.0e9;  %in Pa
E3   = E2; 
G12  = 7.0e9;   %in Pa
G13  = G12;
nu12 = 0.3;
nu13 = nu12;
nu23 = 0.4;

G23 = E2*0.5/(1+nu23);%in Pa

%Layup: [30/45/60]

%-------------------------------------------------------------------
% Step1: Find [S] for material in the material CSYS
[S] = findS(E1, E2, nu12, nu23, G12, G23); % in Pa^(-1)

%-------------------------------------------------------------------
% Step2: Find [S_dash] for each layer in the problem CSYS
Sdash_1 = findSdash(S,30); % in Pa^(-1)
Sdash_2 = findSdash(S,45); % in Pa^(-1)
Sdash_3 = findSdash(S,60); % in Pa^(-1)

%-------------------------------------------------------------------
% Step3: Find Set and Q for each layer
[Se_1, Set_1, St_1] = decomposeS(Sdash_1);
[Se_2, Set_2, St_2] = decomposeS(Sdash_2);
[Se_3, Set_3, St_3] = decomposeS(Sdash_3);

Q_1 = inv(Se_1); Q_2 = inv(Se_2); Q_3 = inv(Se_3);


K1 = (Set_1)'*Q_1; K2 = (Set_2)'*Q_2; K3 = (Set_3)'*Q_3; 

%-------------------------------------------------------------------
% Step4: Formulate fluctuating functions vector for each layer

% --> Step 4.0: Define constants epsilon and kappa
syms e11 e22 e12 real
syms k11 k22 k12 real

e = [e11, e22, 2*e12]';
k = [k11, k22, 2*k12]';

% --> Step 4.1: Make vector of integration constants for each layer
syms c11 c12 c13 real % cij => ith layer, jth constant 
syms c21 c22 c23 real
syms c31 c32 c33 real

c1_vec = [c11; c12; c13];
c2_vec = [c21; c22; c23];
c3_vec = [c31; c32; c33];

% --> Step 4.3: Make the function in x3, replace x3 with just 'X'
syms X
f=  X*e + 0.5*(X^2)*k;

% --> Step 4.4: Find the fluctuating functions for each layer
w_1 = K1*f + c1_vec;
w_2 = K2*f + c2_vec;
w_3 = K3*f + c3_vec;

%% -------------------------------------------------------------------
%Step5: Apply continuity conditions at the interface of each layer

syms t real
z0 = -1.5*t; z1 = -0.5*t; z2 = 0.5*t; z3 = 1.5*t;

% --> Interface 1: Layer 1 and Layer 2
eq1 = subs(w_1,X,z1)-subs(w_2,X,z1)==0;

sols1 = solve(eq1,c21,c22,c23);

% Solve for c_2 in terms of c_1
c21_sol = sols1.c21; c22_sol = sols1.c22; c23_sol = sols1.c23;

% --> Interface 2: Layer 2 and Layer 3
eq2 = subs(subs(w_2,X,z2),{c21,c22,c23},{c21_sol,c22_sol,c23_sol})-subs(w_3,X,z2)==0;

sols2 = solve(eq2,c31,c32,c33);

% Solve for c_2 in terms of c_1
c31_sol = sols2.c31; c32_sol = sols2.c32; c33_sol = sols2.c33;

%-------------------------------------------------------------------
%Step6: Update constants in terms of c_1 for all the layers
w_1_sol = w_1;
w_2_sol = subs(w_2,{c21,c22,c23},{c21_sol,c22_sol,c23_sol});
w_3_sol = subs(w_3,{c31,c32,c33},{c31_sol,c32_sol,c33_sol});


%-------------------------------------------------------------------
%Step7: Average vanishes for w(i)

w_1 = w_1_sol; w_2 = w_2_sol; w_3 = w_3_sol;

w_11 = w_1(1);  w_21 = w_2(1);  w_31 = w_3(1);
w_12 = w_1(2);  w_22 = w_2(2);  w_32 = w_3(2);
w_13 = w_1(3);  w_23 = w_2(3);  w_33 = w_3(3);

% Component w(1) --> average (w_11,w_21,w_31)
avg1 = (1/(3*t))*(int(w_11,X,z0,z1) + int(w_21,X,z1,z2) + int(w_31,X,z2,z3)) == 0;

% Component w(2) --> average (w_12,w_22,w_32)
avg2 = (1/(3*t))*(int(w_12,X,z0,z1) + int(w_22,X,z1,z2) + int(w_32,X,z2,z3)) == 0;

% Component w(3) --> average (w_13,w_23,w_33)
avg3 = (1/(3*t))*(int(w_13,X,z0,z1) + int(w_23,X,z1,z2) + int(w_33,X,z2,z3)) == 0;

%-------------------------------------------------------------------
%Step8: Solve for constants of Layer 1, c11, c12, c13
sols_c1 = solve([avg1,avg2,avg3],[c11,c12,c13]);



%-------------------------------------------------------------------
% Step9: Final values of constants of Layer 1 in terms of t,e,k only
c11_final = sols_c1.c11;
c12_final = sols_c1.c12;
c13_final = sols_c1.c13;

c1_final = [c11_final, c12_final, c13_final]';
% --> Step 9.1: Update constants of all the layers

% Layer 2::
c21_final = subs(c21_sol,{c11,c12,c13},{c11_final,c12_final,c13_final});
c22_final = subs(c22_sol,{c11,c12,c13},{c11_final,c12_final,c13_final});
c23_final = subs(c23_sol,{c11,c12,c13},{c11_final,c12_final,c13_final});

c2_final = [c21_final, c22_final, c23_final]';

% Layer 3::
c31_final = subs(c31_sol,{c11,c12,c13},{c11_final,c12_final,c13_final});
c32_final = subs(c32_sol,{c11,c12,c13},{c11_final,c12_final,c13_final});
c33_final = subs(c33_sol,{c11,c12,c13},{c11_final,c12_final,c13_final});

c3_final = [c31_final, c32_final, c33_final]';

%-------------------------------------------------------------------
% Step 10: Final fluctuating functions in terms of t, e, k

w_1_final = K1*f + c1_final;
w_2_final = K2*f + c2_final;
w_3_final = K3*f + c3_final;




























