function idx=checkFailureNew(f_allow,LAYUP)

[failN22, idx] = min(f_allow);

if failN22 ~= Inf
if idx==1 && LAYUP(1)==abs(LAYUP(2))
    if abs(LAYUP(2))~=abs(LAYUP(3))
        fprintf("Lamina-1 (%0.2f) fails at N22 = %0.2f N/mm \n",LAYUP(1),failN22*1e-3)
        fprintf("Lamina-2 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(2),failN22*1e-3)
    else
        fprintf("Lamina-1 (%0.2f) fails at N22 = %0.2f N/mm \n",LAYUP(1),failN22*1e-3)
        fprintf("Lamina-2 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(2),failN22*1e-3)
        fprintf("Lamina-3 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(3),failN22*1e-3)
        fprintf("Lamina-4 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(4),failN22*1e-3)
    end
elseif idx==1 && LAYUP(1)~=abs(LAYUP(2))
    fprintf("Lamina-1 (%0.2f) fails at N22 = %0.2f N/mm \n",LAYUP(1),failN22*1e-3)
elseif idx==2
    fprintf("Lamina-2 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(2),failN22*1e-3)
elseif idx==3
    if abs(LAYUP(3))==abs(LAYUP(4))
        fprintf("Lamina-3 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(3),failN22*1e-3)
        fprintf("Lamina-4 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(4),failN22*1e-3)
    else
        fprintf("Lamina-3 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(3),failN22*1e-3)
    end
elseif idx==4
    fprintf("Lamina-4 (%0.2f) fails at N22 = %0.2f N/mm\n",LAYUP(4),failN22*1e-3)
else
    fprintf("Error! \n")
end

else
    fprintf("Laminate has already failed \n")
end


end