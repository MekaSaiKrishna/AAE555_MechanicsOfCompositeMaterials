% For Orthotropic Material

function [S] = findS(E1, E2, nu12, nu23, G12, G23)

S = [1/E1, -nu12/E1, -nu12/E1, 0, 0, 0;
    -nu12/E1, 1/E2, -nu23/E2, 0, 0, 0;
    -nu12/E1, -nu23/E2, 1/E2, 0, 0, 0;
    0, 0, 0, 1/G23, 0, 0;
    0, 0, 0, 0, 1/G12, 0;
    0, 0, 0, 0, 0, 1/G12];

end
