
function [E1_star,E2_star,nu12_star,G12_star,eta121_star,eta122_star] = findEffInPlane(Se_star)

Se_star_11  = Se_star(1,1);
Se_star_12  = Se_star(1,2);
Se_star_16  = Se_star(1,3);
Se_star_22  = Se_star(2,2);
Se_star_26  = Se_star(2,3);
Se_star_66  = Se_star(3,3);

E1_star     = 1/(Se_star_11);
E2_star     = 1/(Se_star_22);
nu12_star   = -(Se_star_12)/(Se_star_11);
G12_star    = 1/(Se_star_66);
eta121_star = (Se_star_16)/(Se_star_66);
eta122_star = (Se_star_26)/(Se_star_66);

end