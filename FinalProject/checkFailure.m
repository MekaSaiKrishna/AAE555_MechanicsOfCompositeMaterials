function checkFailure(f_allow)

[failN22, idx] = min(f_allow);

if idx==1
    fprintf("Lamina-1 (+45) fails at N22 = %0.2f N/mm, \n",failN22*1e-3)
% elseif idx==2
    fprintf("Lamina-2 (-45) fails at N22 = %0.2f  N/mm \n",failN22*1e-3)
elseif idx==3
    fprintf("Lamina-3 (0) fails at N22 = %0.2f  N/mm \n",failN22*1e-3)
elseif idx==4
    fprintf("Lamina-4 (90) fails at N22 = %0.2f  N/mm \n",failN22*1e-3)
else
    fprintf("Error! \n")
end

end