function Zs = cons_z_update(Xs, Ys, Zs, rows, row2cp, cpIdx)                

for row = rows
    rIdx = row2cp(row); kIdx = row2cp(cpIdx{row});
    Zs{rIdx} = 0;
    for k = kIdx % neighbors                                           
        Zs{rIdx} = Zs{rIdx} + (Xs{k}(rIdx)+Ys{k}{rIdx})/length(kIdx);
    end
end

end