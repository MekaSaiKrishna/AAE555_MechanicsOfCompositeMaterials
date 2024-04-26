% Project Part-3:

% Layup-1: [+45/-45/0/90]s
% Layup-2: [0/30/60/90]s
% Layup-3: [+45/-45/-45/+45]s
% Layup-4: [0/0/90/90]s

%---------------------------------------------------
% Load HROM engineering constants of lamina for 
% varying volume fiber fraction

% Format: [E1,nu12,E2,nu23,G12]' 
load("consts.mat")



%---------------------------------------------------
% Save Strength Parameters
% Format: [X; Xdash; Y; Ydash; S; R] % in MPa
load("strength.mat")

%---------------------------------------------------
% Save Layup Sequences
LAYUP = [45,-45,  0, 90;  % Layup-1
          0, 30, 60, 90;  % Layup-2
         45,-45,-45, 45;  % Layup-3
          0,  0, 90, 90]; % Layup-4

%---------------------------------------------------

% Design - 1A

% CONSTS and STRENGTH - depend on volume fraction i.e. 'A'
% LAYUP - depends on '1'

vfCode = ["A","B","C"];
vfVals = [0.5,0.6,0.7];

for i=1:1 % Layup Scheme (1,2,3,4)
    for j=1:1 % Volume Fraction (A,B,C)
        fprintf("Design: %.15g %s \n",i,vfCode(j))
        fprintf('Layup: [%s]_s \n', join(string(LAYUP(i,:)), ','));
        fprintf('Volume Fraction: %0.2f \n \n',vfVals(j));
        doMyJob(CONSTS(:,j),LAYUP(i,:),STRENGTH(:,j))
        fprintf("===================================================== \n")
    end
end





















