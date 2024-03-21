function [a,b,d] = find_abd(A,B,D)

    a = inv(A - B*inv(D)*(B'));
    b = -(inv(D)*B*a)';
    d = inv(D - (B')*inv(A)*B);

end
