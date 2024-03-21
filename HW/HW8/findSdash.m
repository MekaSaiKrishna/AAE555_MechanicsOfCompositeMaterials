function Sdash = findSdash(S,t)

Sdash = zeros(6,6);

c = cosd(t); s = sind(t);

S11 = S(1,1);
S12 = S(1,2);
S13 = S(1,3);
S23 = S(2,3);
S22 = S(2,2);
S33 = S(3,3);
S44 = S(4,4);
S55 = S(5,5);
S66 = S(6,6);

Sdash(1,1) = (c^4)*S11 + (c^2*s^2)*(2*S12 + S66) + (s^4)*S22;
Sdash(1,2) = (c^2*s^2)*(S11 + S22 - S66) + (c^4 + s^4)*S12;
Sdash(2,1) = Sdash(1,2);
Sdash(1,3) = c^2*S13 + s^2*S23;
Sdash(3,1) = Sdash(1,3);
Sdash(1,6) = s*c*(c^2*(2*S11 - 2*S12 - S66) + s^2*(2*S12 - 2*S22 + S66));
Sdash(6,1) = Sdash(1,6);

Sdash(2,2) = s^4*S11 + (c^2*s^2)*(2*S12 + S66) + c^4*S22;
Sdash(2,3) = s^2*S13 + c^2*S23;
Sdash(3,2) = Sdash(2,3);
Sdash(2,6) = s*c*(s^2*(2*S11 - 2*S12 - S66) + c^2*(2*S12 - 2*S22 + S66));
Sdash(6,2) = Sdash(2,6);

Sdash(3,3) = S33;
Sdash(3,6) = 2*s*c*(S13-S23);
Sdash(6,3) = Sdash(3,6);

Sdash(4,4) = c^2*S44 + s^2*S55;

Sdash(4,5) = s*c*(S55-S44);
Sdash(5,4) = Sdash(4,5);

Sdash(5,5) = s^2*S44 + c^2*S55;

Sdash(6,6) = 4*c^2*s^2*(S11 - 2*S12 + S22) + ((s^2 - c^2)^2)*S66;




end
