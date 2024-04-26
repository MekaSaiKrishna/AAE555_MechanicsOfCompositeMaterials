function A = findA(Q1,Q2,Q3,Q4)

% For a symmetric laminate with 8 layers

% Step-3: Find [A] matrix
t    = 0.127*1e-3; %in mm (thickness of ply)

syms x real
z0 = -4*t; z1 = -3*t; z2 = -2*t; z3 = -t;         % in mm
z4 = 0;    z5 = t;  z6 = 2*t; z7 = 3*t; z8 = 4*t; % in mm

% Q1 = Qd_p45;  Q2 = Qd_m45;  Q3 = Qd_0;   Q4 = Qd_90;  % in Pa

[A] = vpa(int(Q1,x,z0,z1),4) + vpa(int(Q2,x,z1,z2),4) + vpa(int(Q3,x,z2,z3),4)...
    + vpa(int(Q4,x,z3,z4),4) + vpa(int(Q4,x,z4,z5),4) + ...
      vpa(int(Q3,x,z5,z6),4) + vpa(int(Q2,x,z6,z7),4) + vpa(int(Q1,x,z7,z8),4); % N/m

A = vpa(A,4); % N/m

end