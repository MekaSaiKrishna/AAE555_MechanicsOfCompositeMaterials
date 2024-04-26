
function doMyJob(CONSTS,LAYUP,STRENGTH)

% INPUT:
% CONSTS   - [E1,nu12,E2,nu23,G12] - Lamina props for varying Vf (SI Units)
% LAYUP    - [theta1,theta2,theta3,theta4]s - Angles in degrees
% STRENGTH - [X, Xdash, Y, Ydash, S, R] % in Pa

% Lamina Constants for Vf=0.5
E1_lam   = CONSTS(1); %in Pa
nu12_lam = CONSTS(2);
E2_lam   = CONSTS(3); %in Pa
nu23_lam = CONSTS(4);
G12_lam  = CONSTS(5); %in Pa

nu21_lam = E2_lam*(nu12_lam/E1_lam);

t    = 0.127*1e-3; %in mm (thickness of ply)

% LAYUP(1,:)

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
Q1_0   = findQbar(Q,LAYUP(1));  % in Pa
Q2_0   = findQbar(Q,LAYUP(2));  % in Pa
Q3_0   = findQbar(Q,LAYUP(3));  % in Pa
Q4_0   = findQbar(Q,LAYUP(4));  % in Pa

Qnull = zeros(3,3); %in Pa

%-------------------------------------------------------------------
% Step-3: Find [A] matrix

Q1 = Q1_0;  Q2 = Q2_0;  Q3 = Q3_0;   Q4 = Q4_0;  % in Pa

A = findA(Q1,Q2,Q3,Q4); %in Pa*m

syms N22 real % tensile
eps = inv(A)*[0,N22,0]';

(eps/N22)*10^10;

%-------------------------------------------------------------------
% Step-4: Find Stresses in each layer in the PROBLEM CSYS

% Each column represents layers 1 to 4
% Each row element is sig11, sig22, sig12
stress_prob = sym(zeros(3,4));
Q = sym(zeros(3,3,4));

Q(:,:,1)=Q1; Q(:,:,2)=Q2; Q(:,:,3)=Q3; Q(:,:,4)=Q4;

for i=1:4
stress_prob(:,i) = Q(:,:,i)*eps; %in Pa
end

%-------------------------------------------------------------------
% Step-5: Convert stresses from PROBLEM CSYS to MATERIAL CSYS
stress_mat = sym(zeros(3,4));
R = sym(zeros(3,3,4));

R(:,:,1) = Rsig(LAYUP(1)); R(:,:,2) = Rsig(LAYUP(2)); 
R(:,:,3) = Rsig(LAYUP(3)); R(:,:,4) = Rsig(LAYUP(4));

for i=1:4
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
    f_allow(i) = TsaiWu(STRENGTH(:),stress_mat_by_N22(:,i));
end

fprintf("First Ply/Plies Failure: \n")
idx1 = checkFailureNew(f_allow,LAYUP);
f_allow;
fprintf("----------------------------------------------------- \n")

%-------------------------------------------------------------------
% Step-7: Find degraded [A] after first ply failure


% Case-1: When idx=1
% Account for [+-45/0/90]s
if idx1==1 && LAYUP(1)==abs(LAYUP(2))
    Q1=Qnull;
    Q2=Qnull;
    % Account for [(+-45)s]s
    if abs(LAYUP(2))==abs(LAYUP(3))
        Q3=Qnull;
        Q4=Qnull;
    else
        Q3=Q3_0;
        Q4=Q4_0;
    end
else
    Q1 = Qnull;
    Q2 = Q2_0;
    Q3 = Q3_0;
    Q4 = Q4_0;
end

% Case-2: When idx=2
if idx1==2
    Q1 = Q1_0;
    Q2 = Qnull;
    Q3 = Q3_0;
    Q4 = Q4_0;
end

% Case-3: When idx=3
if idx1==3 
    % Account for [0_2/90_2]s
    if LAYUP(3)==LAYUP(4)
        Q1 = Q1_0;
        Q2 = Q2_0;
        Q3 = Qnull;
        Q4 = Qnull;
    else
        Q1 = Q1_0;
        Q2 = Q2_0;
        Q3 = Qnull;
        Q4 = Q4_0;
    end

end

% Case-4: When idx=4
if idx1==4
    Q1 = Q1_0;
    Q2 = Q2_0;
    Q3 = Q3_0;
    Q4 = Qnull;
end

%-------------------
% Q1
% Q2
% Q3
% Q4

A = findA(Q1,Q2,Q3,Q4); %in N/m

syms N22 real % tensile
eps = inv(A)*[0,N22,0]';


%-------------------------------------------------------------------
% Step-8: Find degraded [A] after first ply failure
A = findA(Q1,Q2,Q3,Q4); %in N/m

syms N22 real % tensile
eps = inv(A)*[0,N22,0]';

% (eps/N22)*10^10

%-------------------------------------------------------------------
% Step-8: Find Stresses in each layer in the PROBLEM CSYS

% Each column represents layers 1 to 4
% Each row element is sig11, sig22, sig12
stress_prob = sym(zeros(3,4));
Q = sym(zeros(3,3,4));

Q(:,:,1)=Q1; Q(:,:,2)=Q2; Q(:,:,3)=Q3; Q(:,:,4)=Q4;


for i=1:4
stress_prob(:,i) = Q(:,:,i)*eps;
end

%-------------------------------------------------------------------
% Step-9: Convert stresses from PROBLEM CSYS to MATERIAL CSYS
stress_mat = sym(zeros(3,4));
R = sym(zeros(3,3,4));

R(:,:,1) = Rsig(LAYUP(1)); R(:,:,2) = Rsig(LAYUP(2)); 
R(:,:,3) = Rsig(LAYUP(3)); R(:,:,4) = Rsig(LAYUP(4));

for i=1:4
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

fprintf("Second Ply/Plies Failure: \n")
idx1 = checkFailureNew(f_allow,LAYUP);
f_allow;
fprintf("----------------------------------------------------- \n")


%-------------------------------------------------------------------
% Step-11: Find degraded [A] after first ply failure


% Case-1: When idx=1
% Account for [+-45/0/90]s
if idx1==1 && LAYUP(1)==abs(LAYUP(2))
    Q1=Qnull;
    Q2=Qnull;
    % Account for [(+-45)s]s
    if abs(LAYUP(2))==abs(LAYUP(3))
        Q3=Qnull;
        Q4=Qnull;
    %else
    %    Q3=Q3_0;
    %    Q4=Q4_0;
    end
else
    Q1 = Qnull;
    %Q2 = Q2_0;
    %Q3 = Q3_0;
    %Q4 = Q4_0;
end

% Case-2: When idx=2
if idx1==2
    %Q1 = Q1_0;
    Q2 = Qnull;
    %Q3 = Q3_0;
    %Q4 = Q4_0;
end

% Case-3: When idx=3
if idx1==3 
    % Account for [0_2/90_2]s
    if LAYUP(3)==LAYUP(4)
        %Q1 = Q1_0;
        %Q2 = Q2_0;
        Q3 = Qnull;
        Q4 = Qnull;
    else
        %Q1 = Q1_0;
        %Q2 = Q2_0;
        Q3 = Qnull;
        %Q4 = Q4_0;
    end

end

% Case-4: When idx=4
if idx1==4
    %Q1 = Q1_0;
    %Q2 = Q2_0;
    %Q3 = Q3_0;
    Q4 = Qnull;
end

%-------------------
% Q1
% Q2
% Q3
% Q4

A = findA(Q1,Q2,Q3,Q4); %in N/m

syms N22 real % tensile
eps = inv(A)*[0,N22,0]';


%-------------------------------------------------------------------
% Step-12: Find Stresses in each layer in the PROBLEM CSYS

% Each column represents layers 1 to 4
% Each row element is sig11, sig22, sig12
stress_prob = sym(zeros(3,4));
Q = sym(zeros(3,3,4));

Q(:,:,1)=Q1; Q(:,:,2)=Q2; Q(:,:,3)=Q3; Q(:,:,4)=Q4;


for i=1:4
stress_prob(:,i) = Q(:,:,i)*eps;
end

%-------------------------------------------------------------------
% Step-13: Convert stresses from PROBLEM CSYS to MATERIAL CSYS
stress_mat = sym(zeros(3,4));
R = sym(zeros(3,3,4));

R(:,:,1) = Rsig(LAYUP(1)); R(:,:,2) = Rsig(LAYUP(2)); 
R(:,:,3) = Rsig(LAYUP(3)); R(:,:,4) = Rsig(LAYUP(4));

for i=1:4
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

fprintf("Third Ply/Plies Failure: \n")
idx1 = checkFailureNew(f_allow,LAYUP);
f_allow;
fprintf("----------------------------------------------------- \n")


%-------------------------------------------------------------------
% Step-15: Find degraded [A] after first ply failure


% Case-1: When idx=1
% Account for [+-45/0/90]s
if idx1==1 && LAYUP(1)==abs(LAYUP(2))
    Q1=Qnull;
    Q2=Qnull;
    % Account for [(+-45)s]s
    if abs(LAYUP(2))==abs(LAYUP(3))
        Q3=Qnull;
        Q4=Qnull;
    %else
    %    Q3=Q3_0;
    %    Q4=Q4_0;
    end
else
    Q1 = Qnull;
    %Q2 = Q2_0;
    %Q3 = Q3_0;
    %Q4 = Q4_0;
end

% Case-2: When idx=2
if idx1==2
    %Q1 = Q1_0;
    Q2 = Qnull;
    %Q3 = Q3_0;
    %Q4 = Q4_0;
end

% Case-3: When idx=3
if idx1==3 
    % Account for [0_2/90_2]s
    if LAYUP(3)==LAYUP(4)
        %Q1 = Q1_0;
        %Q2 = Q2_0;
        Q3 = Qnull;
        Q4 = Qnull;
    else
        %Q1 = Q1_0;
        %Q2 = Q2_0;
        Q3 = Qnull;
        %Q4 = Q4_0;
    end

end

% Case-4: When idx=4
if idx1==4
    %Q1 = Q1_0;
    %Q2 = Q2_0;
    %Q3 = Q3_0;
    Q4 = Qnull;
end

%-------------------
% Q1
% Q2
% Q3
% Q4

A = findA(Q1,Q2,Q3,Q4); %in N/m

syms N22 real % tensile
eps = inv(A)*[0,N22,0]';


%-------------------------------------------------------------------
% Step-16: Find Stresses in each layer in the PROBLEM CSYS

% Each column represents layers 1 to 4
% Each row element is sig11, sig22, sig12
stress_prob = sym(zeros(3,4));
Q = sym(zeros(3,3,4));

Q(:,:,1)=Q1; Q(:,:,2)=Q2; Q(:,:,3)=Q3; Q(:,:,4)=Q4;


for i=1:4
stress_prob(:,i) = Q(:,:,i)*eps;
end

%-------------------------------------------------------------------
% Step-17: Convert stresses from PROBLEM CSYS to MATERIAL CSYS
stress_mat = sym(zeros(3,4));
R = sym(zeros(3,3,4));

R(:,:,1) = Rsig(LAYUP(1)); R(:,:,2) = Rsig(LAYUP(2)); 
R(:,:,3) = Rsig(LAYUP(3)); R(:,:,4) = Rsig(LAYUP(4));

for i=1:4
stress_mat(:,i) = inv(R(:,:,i))*stress_prob(:,i); %in Pa
end

% Stress in material coordinate system divided by N22
stress_mat_by_N22 = vpa(stress_mat/N22, 5);

%-------------------------------------------------------------------
% Step-18: Apply Tsai-Wu Failure Criteria for each layer

% Store allowable N22 for each layer
f_allow=zeros(1,4); 

for i=1:4
    % For lamina with Vf=0.5
    f_allow(i) = TsaiWu(STRENGTH(:,1),stress_mat_by_N22(:,i));
end

fprintf("Fourth Ply/Plies Failure: \n")
idx1 = checkFailureNew(f_allow,LAYUP);
f_allow;


end