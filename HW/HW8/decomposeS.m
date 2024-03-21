function [Se, Set, St] = decomposeS(S)

Se = [S(1,1),S(1,2),S(1,6);
      S(1,2),S(2,2),S(2,6);
      S(1,6),S(2,6),S(6,6)];

Set = [S(1,3),S(1,4),S(1,5);
       S(2,3),S(2,4),S(2,5);
       S(3,6),S(4,6),S(5,6)];

St = [S(3,3),S(3,4),S(3,5);
      S(3,4),S(4,4),S(4,5);
      S(3,5),S(4,5),S(5,5)];

end
