% Homework -10: Question 5.7

% Layup: [+45/-45/0]s => [+45/-45/0/0/-45/+45]

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

t = 0.127*1e-3; %in m (thickness of ply)

%%
% (i) Compute the effective in-plane properties 
% in terms of engineering constants using Eq. (5.76) and compare 
% the results with those computed using the hybrid rules of mixtures 
% for laminates in Section 4.6.2

% Find [A,B,D] matrix

% Find: Se_star = h*inv(A)

% Find Effec In-plane props from Se_star

%------------------------------------------------
% USE HROM MATLAB CODE FROM PREVIOUS HOMEWORK

% Step-1: Find Q for each Layer

% Step-2: Find Qstar 

% Step-3: Se_star = inv(Qstar)

% Step-4: Compute effective engineering constants from Se_star 



% Note: assume any value for the rest of the material properties needed 
% for the hybrid rules of mixtures and the effective in-plane properties 
% will not be affected by this arbitrariness.

























