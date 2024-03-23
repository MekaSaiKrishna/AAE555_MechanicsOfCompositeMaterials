% Question 5.5

% Layup: [0/60/-60]s => [0/60/-60/-60/60/0]

% Input Properties:
E1   = 140e9; %in Pa = N/m^2
E2   = 10e9;  %in Pa = N/m^2
G12  = 7e9;   %in Pa = N/m^2
nu12 = 0.3;

nu21 = E2*(nu12/E1);

t    = 0.127*1e-3; %in m

%-------------------------------------------------------------------
% Step1: Find [Q] for material in the material CSYS
Q11 = E1/(1-nu12*nu21);
Q22 = E2/(1-nu12*nu21);
Q12 = nu12*E2/(1-nu12*nu21);
Q66 = G12;

Q = [Q11, Q12, 0;
    Q12, Q22, 0;
    0, 0, Q66] %in Pa

%-------------------------------------------------------------------
% Step-2: Find Qbar for different angles
[Qd_0]   = findQbar(Q,0)  %in Pa = N/m^2
[Qd_p60] = findQbar(Q,60) %in Pa = N/m^2
[Qd_m60] = findQbar(Q,-60) %in Pa = N/m^2


%-------------------------------------------------------------------
% Step-3: Find [A, B, D] matrix
syms x real
z0 = -3*t; z1 = -2*t; z2 = -t; z3 = 0; %in m
z4 = t;    z5 = 2*t;  z6 = 3*t; %in m 

%Layup: [0/60/-60]s => [0/60/-60/-60/60/0]

Q1 = Qd_0;   Q2 = Qd_p60;  Q3 = Qd_m60; %in Pa
Q4 = Qd_m60; Q5 = Qd_p60;  Q6 = Qd_0;   %in Pa


[A] = vpa(int(Q1,x,z0,z1),4) + vpa(int(Q2,x,z1,z2),4) + vpa(int(Q3,x,z2,z3),4)...
    + vpa(int(Q4,x,z3,z4),4) + vpa(int(Q5,x,z4,z5),4) + ...
      vpa(int(Q6,x,z5,z6),4); %in = (N/m^2)*m = (N/m)

A = vpa(A,4) %in = (N/m^2)*m = (N/m)

[B] = vpa(int(Q1*x,x,z0,z1),4) + vpa(int(Q2*x,x,z1,z2),4) + vpa(int(Q3*x,x,z2,z3),4)...
    + vpa(int(Q4*x,x,z3,z4),4) + vpa(int(Q5*x,x,z4,z5),4) + ...
      vpa(int(Q6*x,x,z5,z6),4); %in Pa*m^2 = (N/m^2)*(m)^2  = N

B = vpa(B,4) %in Pa*m^2 = (N/m^2)*(m)^2  = N


[D] = vpa(int(Q1*x^2,x,z0,z1),4) + vpa(int(Q2*x^2,x,z1,z2),4) + vpa(int(Q3*x^2,x,z2,z3),4)...
    + vpa(int(Q4*x^2,x,z3,z4),4) + vpa(int(Q5*x^2,x,z4,z5),4) + ...
      vpa(int(Q6*x^2,x,z5,z6),4); %in Pa*m^3 = (N/m^2)*(m)^3 = (N*m)

D = vpa(D,4) %in Pa*m^3 = (N/m^2)*(m)^3 = (N*m)

%-------------------------------------------------------------------
% Step-4: Find [a, b, d] matrix
[a,b,d] = find_abd(A,B,D);

abd = [a, b;
       b',d]; % SI Units

%-------------------------------------------------------------------
% Step-5: Loading condition

%==================================================================
% Loads:

% Load-1:
N22 = 1e4; % in N/m

%NM_vec = [N11, N22, N12, M11, M22, M12];
NM_vec1 = [0, N22, 0, 0, 0, 0]';

%ek_vec = [E11, E22, 2*E12, k11, k22, 2*k12];


%------------------------------------------------------------------
% Load-2:
M12 = 10;  % in N

%NM_vec = [N11, N22, N12, M11, M22, M12];
NM_vec2 = [0, 0, 0, 0, 0, M12]';

%ek_vec = [E11, E22, 2*E12, k11, k22, 2*k12];

%-------------------------------------------------------------------
% Step-6: Plate Constitutive Relation

% Load-1:
ek_vec1 = abd*NM_vec1;
e_vec1 = ek_vec1(1:3) 
k_vec1 = ek_vec1(4:6)

syms x real
e_load1 = e_vec1 + x*k_vec1;

% Layer 1: z0 <= x <= z1
s_load1_1 = Q1*e_load1;
x1 = linspace(z0,z1,10);

% Layer 2: z1 <= x <= z2
s_load1_2 = Q2*e_load1;
x2 = linspace(z1,z2,10);

% Layer 3: z2 <= x <= z3
s_load1_3 = Q3*e_load1;
x3 = linspace(z2,z3,10);

% Layer 4: z3 <= x <= z4
s_load1_4 = Q4*e_load1;
x4 = linspace(z3,z4,10);

% Layer 5: z4 <= x <= z5
s_load1_5 = Q5*e_load1;
x5 = linspace(z4,z5,10);

% Layer 6: z5 <= x <= z6
s_load1_6 = Q6*e_load1;
x6 = linspace(z5,z6,10);

% 1) Plot sigma_11

figure()

%Layer-1
plot(subs(s_load1_1(1),x,x1)*1e-6,x1*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load1_1(1),x,x1(end)),subs(s_load1_2(1),x,x2(1)) ]*1e-6,[x1(end),x2(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-2
plot(subs(s_load1_2(1),x,x2)*1e-6,x2*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load1_2(1),x,x2(end)),subs(s_load1_3(1),x,x3(1)) ]*1e-6,[x2(end),x3(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-3
plot(subs(s_load1_3(1),x,x3)*1e-6,x3*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load1_3(1),x,x3(end)),subs(s_load1_4(1),x,x4(1)) ]*1e-6,[x3(end),x4(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-4
plot(subs(s_load1_4(1),x,x4)*1e-6,x4*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load1_4(1),x,x4(end)),subs(s_load1_5(1),x,x5(1)) ]*1e-6,[x4(end),x5(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-5
plot(subs(s_load1_5(1),x,x5)*1e-6,x5*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load1_5(1),x,x5(end)),subs(s_load1_6(1),x,x6(1)) ]*1e-6,[x5(end),x6(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-6
plot(subs(s_load1_6(1),x,x6)*1e-6,x6*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold off
[t,s] = title('\sigma_{11} vs z','Load 1: N22 = 1e4 N/m',...
    'Color','blue');
xlabel("\sigma_{11} (in MPa)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])

%----------------------------------------------------------------------------------------------
% 2) Plot sigma_22
figure()

%Layer-1
plot(subs(s_load1_1(2),x,x1)*1e-6,x1*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load1_1(2),x,x1(end)),subs(s_load1_2(2),x,x2(1)) ]*1e-6,[x1(end),x2(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-2
plot(subs(s_load1_2(2),x,x2)*1e-6,x2*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load1_2(2),x,x2(end)),subs(s_load1_3(2),x,x3(1)) ]*1e-6,[x2(end),x3(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-3
plot(subs(s_load1_3(2),x,x3)*1e-6,x3*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load1_3(2),x,x3(end)),subs(s_load1_4(2),x,x4(1)) ]*1e-6,[x3(end),x4(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-4
plot(subs(s_load1_4(2),x,x4)*1e-6,x4*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load1_4(2),x,x4(end)),subs(s_load1_5(2),x,x5(1)) ]*1e-6,[x4(end),x5(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-5
plot(subs(s_load1_5(2),x,x5)*1e-6,x5*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load1_5(2),x,x5(end)),subs(s_load1_6(2),x,x6(1)) ]*1e-6,[x5(end),x6(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-6
plot(subs(s_load1_6(2),x,x6)*1e-6,x6*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold off
[t,s] = title('\sigma_{22} vs z','Load 1: N22 = 1e4 N/m',...
    'Color','blue');
xlabel("\sigma_{22} (in MPa)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])

%----------------------------------------------------------------------------------------------
% 3) Plot sigma_12
figure()

%Layer-1
plot(subs(s_load1_1(3),x,x1)*1e-6,x1*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load1_1(3),x,x1(end)),subs(s_load1_2(3),x,x2(1)) ]*1e-6,[x1(end),x2(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-2
plot(subs(s_load1_2(3),x,x2)*1e-6,x2*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load1_2(3),x,x2(end)),subs(s_load1_3(3),x,x3(1)) ]*1e-6,[x2(end),x3(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-3
plot(subs(s_load1_3(3),x,x3)*1e-6,x3*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load1_3(3),x,x3(end)),subs(s_load1_4(3),x,x4(1)) ]*1e-6,[x3(end),x4(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-4
plot(subs(s_load1_4(3),x,x4)*1e-6,x4*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load1_4(3),x,x4(end)),subs(s_load1_5(3),x,x5(1)) ]*1e-6,[x4(end),x5(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-5
plot(subs(s_load1_5(3),x,x5)*1e-6,x5*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load1_5(3),x,x5(end)),subs(s_load1_6(3),x,x6(1)) ]*1e-6,[x5(end),x6(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-6
plot(subs(s_load1_6(3),x,x6)*1e-6,x6*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold off
[t,s] = title('\sigma_{12} vs z','Load 1: N22 = 1e4 N/m',...
    'Color','blue');
xlabel("\sigma_{12} (in MPa)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])


%------------------------------------------------------------------------
% 4) Plot epsilon_11

figure()

%Layer-1
plot(subs(e_load1(1),x,x1)*1e6,x1*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-2
plot(subs(e_load1(1),x,x2)*1e6,x2*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-3
plot(subs(e_load1(1),x,x3)*1e6,x3*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-4
plot(subs(e_load1(1),x,x4)*1e6,x4*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-5
plot(subs(e_load1(1),x,x5)*1e6,x5*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-6
plot(subs(e_load1(1),x,x6)*1e6,x6*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold off
[t,s] = title('\epsilon_{11} vs z','Load 1: N22 = 1e4 N/m',...
    'Color','blue');
xlabel("\epsilon_{11} (micro strains)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])

%------------------------------------------------------------------------
% 5) Plot epsilon_22

figure()

%Layer-1
plot(subs(e_load1(2),x,x1)*1e6,x1*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-2
plot(subs(e_load1(2),x,x2)*1e6,x2*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-3
plot(subs(e_load1(2),x,x3)*1e6,x3*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-4
plot(subs(e_load1(2),x,x4)*1e6,x4*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-5
plot(subs(e_load1(2),x,x5)*1e6,x5*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-6
plot(subs(e_load1(2),x,x6)*1e6,x6*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold off
[t,s] = title('\epsilon_{22} vs z','Load 1: N22 = 1e4 N/m',...
    'Color','blue');
xlabel("\epsilon_{22} (micro strains)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])

%------------------------------------------------------------------------
% 6) Plot epsilon_12

figure()

%Layer-1
plot(subs(0.5*e_load1(3),x,x1)*1e6,x1*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-2
plot(subs(0.5*e_load1(3),x,x2)*1e6,x2*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-3
plot(subs(0.5*e_load1(3),x,x3)*1e6,x3*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-4
plot(subs(0.5*e_load1(3),x,x4)*1e6,x4*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-5
plot(subs(0.5*e_load1(3),x,x5)*1e6,x5*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-6
plot(subs(0.5*e_load1(3),x,x6)*1e6,x6*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold off
[t,s] = title('\epsilon_{12} vs z','Load 1: N22 = 1e4 N/m',...
    'Color','blue');
xlabel("\epsilon_{12} (micro strains)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])


%================================================
% Load-2:
ek_vec2 = abd*NM_vec2;
e_vec2 = ek_vec2(1:3)
k_vec2 = ek_vec2(4:6)

syms x real
e_load2 = e_vec2 + x*k_vec2;

% Layer 1: z0 <= x <= z1
s_load2_1 = Q1*e_load2;

% Layer 2: z1 <= x <= z2
s_load2_2 = Q2*e_load2;

% Layer 3: z2 <= x <= z3
s_load2_3 = Q3*e_load2;

% Layer 4: z3 <= x <= z4
s_load2_4 = Q4*e_load2;

% Layer 5: z4 <= x <= z5
s_load2_5 = Q5*e_load2;

% Layer 6: z5 <= x <= z6
s_load2_6 = Q6*e_load2;


% 1) Plot sigma_11

figure()

%Layer-1
plot(subs(s_load2_1(1),x,x1)*1e-6,x1*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load2_1(1),x,x1(end)),subs(s_load2_2(1),x,x2(1)) ]*1e-6,[x1(end),x2(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-2
plot(subs(s_load2_2(1),x,x2)*1e-6,x2*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load2_2(1),x,x2(end)),subs(s_load2_3(1),x,x3(1)) ]*1e-6,[x2(end),x3(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-3
plot(subs(s_load2_3(1),x,x3)*1e-6,x3*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load2_3(1),x,x3(end)),subs(s_load2_4(1),x,x4(1)) ]*1e-6,[x3(end),x4(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-4
plot(subs(s_load2_4(1),x,x4)*1e-6,x4*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load2_4(1),x,x4(end)),subs(s_load2_5(1),x,x5(1)) ]*1e-6,[x4(end),x5(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-5
plot(subs(s_load2_5(1),x,x5)*1e-6,x5*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

plot([subs(s_load2_5(1),x,x5(end)),subs(s_load2_6(1),x,x6(1)) ]*1e-6,[x5(end),x6(1)]*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-6
plot(subs(s_load2_6(1),x,x6)*1e-6,x6*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold off
[t,s] = title('\sigma_{11} vs z','Load-2: M12 = 10 N',...
    'Color','blue');
xlabel("\sigma_{11} (in MPa)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])


%----------------------------------------------------------------------------------------------
% 2) Plot sigma_22
figure()

%Layer-1
plot(subs(s_load2_1(2),x,x1)*1e-6,x1*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load2_1(2),x,x1(end)),subs(s_load2_2(2),x,x2(1)) ]*1e-6,[x1(end),x2(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-2
plot(subs(s_load2_2(2),x,x2)*1e-6,x2*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load2_2(2),x,x2(end)),subs(s_load2_3(2),x,x3(1)) ]*1e-6,[x2(end),x3(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-3
plot(subs(s_load2_3(2),x,x3)*1e-6,x3*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load2_3(2),x,x3(end)),subs(s_load2_4(2),x,x4(1)) ]*1e-6,[x3(end),x4(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-4
plot(subs(s_load2_4(2),x,x4)*1e-6,x4*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load2_4(2),x,x4(end)),subs(s_load2_5(2),x,x5(1)) ]*1e-6,[x4(end),x5(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-5
plot(subs(s_load2_5(2),x,x5)*1e-6,x5*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

plot([subs(s_load2_5(2),x,x5(end)),subs(s_load2_6(2),x,x6(1)) ]*1e-6,[x5(end),x6(1)]*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-6
plot(subs(s_load2_6(2),x,x6)*1e-6,x6*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold off
[t,s] = title('\sigma_{22} vs z','Load-2: M12 = 10 N',...
    'Color','blue');
xlabel("\sigma_{22} (in MPa)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])

%----------------------------------------------------------------------------------------------
% 3) Plot sigma_12
figure()

%Layer-1
plot(subs(s_load2_1(3),x,x1)*1e-6,x1*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load2_1(3),x,x1(end)),subs(s_load2_2(3),x,x2(1)) ]*1e-6,[x1(end),x2(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-2
plot(subs(s_load2_2(3),x,x2)*1e-6,x2*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load2_2(3),x,x2(end)),subs(s_load2_3(3),x,x3(1)) ]*1e-6,[x2(end),x3(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-3
plot(subs(s_load2_3(3),x,x3)*1e-6,x3*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load2_3(3),x,x3(end)),subs(s_load2_4(3),x,x4(1)) ]*1e-6,[x3(end),x4(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-4
plot(subs(s_load2_4(3),x,x4)*1e-6,x4*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load2_4(3),x,x4(end)),subs(s_load2_5(3),x,x5(1)) ]*1e-6,[x4(end),x5(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-5
plot(subs(s_load2_5(3),x,x5)*1e-6,x5*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

plot([subs(s_load2_5(3),x,x5(end)),subs(s_load2_6(3),x,x6(1)) ]*1e-6,[x5(end),x6(1)]*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-6
plot(subs(s_load2_6(3),x,x6)*1e-6,x6*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold off
[t,s] = title('\sigma_{12} vs z','Load-2: M12 = 10 N',...
    'Color','blue');
xlabel("\sigma_{12} (in MPa)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])

%------------------------------------------------------------------------
% 4) Plot epsilon_11

figure()

%Layer-1
plot(subs(e_load2(1),x,x1)*1e6,x1*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-2
plot(subs(e_load2(1),x,x2)*1e6,x2*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-3
plot(subs(e_load2(1),x,x3)*1e6,x3*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-4
plot(subs(e_load2(1),x,x4)*1e6,x4*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-5
plot(subs(e_load2(1),x,x5)*1e6,x5*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold on

%Layer-6
plot(subs(e_load2(1),x,x6)*1e6,x6*1e3,'Color',[0 0.4470 0.7410],LineWidth=2.5)
hold off
[t,s] = title('\epsilon_{11} vs z','Load 2: M12 = 10 N',...
    'Color','blue');
xlabel("\epsilon_{11} (micro strains)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])

%------------------------------------------------------------------------
% 5) Plot epsilon_22

figure()

%Layer-1
plot(subs(e_load2(2),x,x1)*1e6,x1*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-2
plot(subs(e_load2(2),x,x2)*1e6,x2*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-3
plot(subs(e_load2(2),x,x3)*1e6,x3*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-4
plot(subs(e_load2(2),x,x4)*1e6,x4*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-5
plot(subs(e_load2(2),x,x5)*1e6,x5*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold on

%Layer-6
plot(subs(e_load2(2),x,x6)*1e6,x6*1e3,'Color',[0.8500 0.3250 0.0980],LineWidth=2.5)
hold off
[t,s] = title('\epsilon_{22} vs z','Load 2: M12 = 10 N',...
    'Color','blue');
xlabel("\epsilon_{22} (micro strains)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])

%------------------------------------------------------------------------
% 6) Plot epsilon_12

figure()

%Layer-1
plot(subs(0.5*e_load2(3),x,x1)*1e6,x1*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-2
plot(subs(0.5*e_load2(3),x,x2)*1e6,x2*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-3
plot(subs(0.5*e_load2(3),x,x3)*1e6,x3*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-4
plot(subs(0.5*e_load2(3),x,x4)*1e6,x4*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-5
plot(subs(0.5*e_load2(3),x,x5)*1e6,x5*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold on

%Layer-6
plot(subs(0.5*e_load2(3),x,x6)*1e6,x6*1e3,'Color',[0.4660 0.6740 0.1880],LineWidth=2.5)
hold off
[t,s] = title('\epsilon_{12} vs z','Load 2: M12 = 10 N',...
    'Color','blue');
xlabel("\epsilon_{12} (micro strains)")
ylabel("z (in mm)")
ylim([z0*1e3,z6*1e3])
