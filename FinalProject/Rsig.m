function R = Rsig(theta)

c = cosd(theta); s = sind(theta);

R    = [c^2,  s^2, -2*s*c;
        s^2,  c^2,  2*s*c;
        s*c, -s*c, c^2-s^2];

end