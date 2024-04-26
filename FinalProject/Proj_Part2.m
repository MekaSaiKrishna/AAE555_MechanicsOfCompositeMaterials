% Project Part-2


% Save Strength Parameters
% Format: [X; Xdash; Y; Ydash; S; R] % in Pa


%     Vf - 0.5, 0.6, 0.7
STRENGTH = [2100, 2400, 2500;
            2100, 2400, 2500;
              62,   68,   71;
             210,  230,  240;
             100,  110,  120;
             100,  110,  120]*1e6; %in Pa

save("strength.mat","STRENGTH")

%---------------------------------------------------
% Load HROM engineering constants of lamina for 
% varying volume fiber fraction
load("consts.mat")

% Lamina Constants for Vf=0.5
E1_lam   = CONSTS(1,1);
nu12_lam = CONSTS(2,1);
E2_lam   = CONSTS(3,1);
nu23_lam = CONSTS(4,1);
G12_lam  = CONSTS(5,1);

nu21_lam = E2_lam*(nu12_lam/E1_lam);

t    = 0.127*1e-3; %in mm (thickness of ply)

%-------------------------------------------------------------------
% Step1: Find [Q] for material in the material CSYS
Q11 = E1_lam/(1-nu12_lam*nu21_lam);
Q22 = E2_lam/(1-nu12_lam*nu21_lam);
Q12 = nu12_lam*E2_lam/(1-nu12_lam*nu21_lam);
Q66 = G12_lam;

Q = [Q11, Q12, 0;
    Q12, Q22, 0;
    0, 0, Q66]; %in Pa

%-------------------------------------------------------------------
% Step-2: Find Qbar for different angles
[Qd_0]   = findQbar(Q,0);   % in Pa*m
[Qd_p45] = findQbar(Q,45);  % in Pa*m
[Qd_m45] = findQbar(Q,-45); % in Pa*m
[Qd_90]  = findQbar(Q,90);  % in Pa*m

%-------------------------------------------------------------------
% Step-3: Find [A] matrix

Q1 = Qd_p45;  Q2 = Qd_m45;  Q3 = Qd_0;   Q4 = Qd_90;  % in Pa

Qnull = zeros(3,3); %in Pa

A = findA(Q1,Q2,Q3,Q4); %in Pa*m

syms N22 real % tensile
eps = inv(A)*[0,N22,0]';

(eps/N22)*10^10;


%-------------------------------------------------------------------
% Step-4: Find Stresses in each layer in the PROBLEM CSYS

% Each column represents layers 1 to 8
% Each row element is sig11, sig22, sig12
stress_prob = sym(zeros(3,8));
Q = sym(zeros(3,3,8));

Q(:,:,1)=Q1; Q(:,:,2)=Q2; Q(:,:,3)=Q3; Q(:,:,4)=Q4;
Q(:,:,5)=Q4; Q(:,:,6)=Q3; Q(:,:,7)=Q2; Q(:,:,8)=Q1;

for i=1:8
stress_prob(:,i) = Q(:,:,i)*eps; %in Pa
end

%-------------------------------------------------------------------
% Step-5: Convert stresses from PROBLEM CSYS to MATERIAL CSYS
stress_mat = sym(zeros(3,8));
R = sym(zeros(3,3,8));

R(:,:,1) = Rsig(+45); R(:,:,2) = Rsig(-45); R(:,:,3) = Rsig(0);   R(:,:,4) = Rsig(90);
R(:,:,5) = Rsig(90);  R(:,:,6) = Rsig(0);   R(:,:,7) = Rsig(-45); R(:,:,8) = Rsig(+45);

for i=1:8
stress_mat(:,i) = inv(R(:,:,i))*stress_prob(:,i); %in Pa
end

% Stress in material coordinate system divided by N22
stress_mat_by_N22 = vpa(stress_mat/N22, 5);

%-------------------------------------------------------------------
% Step-6: Apply Tsai-Wu Failure Criteria for each layer

% Store allowable N22 for each layer
f_allow=zeros(1,4); 

for i=1:4
    % For lamina with Vf=0.5
    f_allow(i) = TsaiWu(STRENGTH(:,1),stress_mat_by_N22(:,i));
end

fprintf("First Ply Failure \n")
checkFailure(f_allow)
f_allow;
fprintf("------------------------------------------- \n")

%% -------------------------------------------------------------------
% Step-7: Find degraded [A] after first ply failure

Q1 = Qd_p45;  Q2 = Qd_m45;  Q3 = Qnull;   Q4 = Qd_90;  % in Pa

A = findA(Q1,Q2,Q3,Q4); %in N/m

syms N22 real % tensile
eps = inv(A)*[0,N22,0]';

(eps/N22)*10^10;

%-------------------------------------------------------------------
% Step-8: Find Stresses in each layer in the PROBLEM CSYS

% Each column represents layers 1 to 8
% Each row element is sig11, sig22, sig12
stress_prob = sym(zeros(3,8));
Q = sym(zeros(3,3,8));

Q(:,:,1)=Q1; Q(:,:,2)=Q2; Q(:,:,3)=Q3; Q(:,:,4)=Q4;
Q(:,:,5)=Q4; Q(:,:,6)=Q3; Q(:,:,7)=Q2; Q(:,:,8)=Q1;

for i=1:8
stress_prob(:,i) = Q(:,:,i)*eps;
end

%-------------------------------------------------------------------
% Step-9: Convert stresses from PROBLEM CSYS to MATERIAL CSYS
stress_mat = sym(zeros(3,8));
R = sym(zeros(3,3,8));

R(:,:,1) = Rsig(+45); R(:,:,2) = Rsig(-45); R(:,:,3) = Rsig(0);   R(:,:,4) = Rsig(90);
R(:,:,5) = Rsig(90);  R(:,:,6) = Rsig(0);   R(:,:,7) = Rsig(-45); R(:,:,8) = Rsig(+45);

for i=1:8
stress_mat(:,i) = inv(R(:,:,i))*stress_prob(:,i); %in Pa
end

% Stress in material coordinate system divided by N22
stress_mat_by_N22 = vpa(stress_mat/N22, 5);

%-------------------------------------------------------------------
% Step-10: Apply Tsai-Wu Failure Criteria for each layer

% Store allowable N22 for each layer
f_allow=zeros(1,4); 

for i=1:4
    % For lamina with Vf=0.5
    f_allow(i) = TsaiWu(STRENGTH(:,1),stress_mat_by_N22(:,i));
end

fprintf("Second Ply/Plies Failure \n")
checkFailure(f_allow)
f_allow;
fprintf("------------------------------------------- \n")

%% -------------------------------------------------------------------
% Step-11: Find degraded [A] after first ply failure

Q1 = Qnull;  Q2 = Qnull;  Q3 = Qnull;   Q4 = Qd_90;  % in Pa

A = findA(Q1,Q2,Q3,Q4); %in N/m

syms N22 real % tensile
eps = inv(A)*[0,N22,0]';

(eps/N22)*10^10;

%-------------------------------------------------------------------
% Step-12: Find Stresses in each layer in the PROBLEM CSYS

% Each column represents layers 1 to 8
% Each row element is sig11, sig22, sig12
stress_prob = sym(zeros(3,8));
Q = sym(zeros(3,3,8));

Q(:,:,1)=Q1; Q(:,:,2)=Q2; Q(:,:,3)=Q3; Q(:,:,4)=Q4;
Q(:,:,5)=Q4; Q(:,:,6)=Q3; Q(:,:,7)=Q2; Q(:,:,8)=Q1;

for i=1:8
stress_prob(:,i) = Q(:,:,i)*eps;
end

%-------------------------------------------------------------------
% Step-13: Convert stresses from PROBLEM CSYS to MATERIAL CSYS
stress_mat = sym(zeros(3,8));
R = sym(zeros(3,3,8));

R(:,:,1) = Rsig(+45); R(:,:,2) = Rsig(-45); R(:,:,3) = Rsig(0);   R(:,:,4) = Rsig(90);
R(:,:,5) = Rsig(90);  R(:,:,6) = Rsig(0);   R(:,:,7) = Rsig(-45); R(:,:,8) = Rsig(+45);

for i=1:8
stress_mat(:,i) = inv(R(:,:,i))*stress_prob(:,i); %in Pa
end

% Stress in material coordinate system divided by N22
stress_mat_by_N22 = vpa(stress_mat/N22, 5);


%-------------------------------------------------------------------
% Step-14: Apply Tsai-Wu Failure Criteria for each layer

% Store allowable N22 for each layer
f_allow=zeros(1,4); 

for i=1:4
    % For lamina with Vf=0.5
    f_allow(i) = TsaiWu(STRENGTH(:,1),stress_mat_by_N22(:,i));
end

fprintf("Third Ply Failure \n")
checkFailure(f_allow)
f_allow;
fprintf("------------------------------------------- \n")





















































