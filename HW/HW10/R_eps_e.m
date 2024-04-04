function R = R_eps_e(t)

c = cosd(t); s = sind(t);

R = [c^2, s^2, -s*c;
    s^2, c^2, s*c;
    2*s*c, -2*s*c, c^2 - s^2];

end