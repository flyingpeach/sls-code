function Zs = cons_z_update(Xs, Ys, Zs, rows, cpIdx)                

for row = rows
    Zs{row} = 0;
    for k = cpIdx{row} % neighbors                                           
        Zs{row} = Zs{row} + (Xs{k}(row)+Ys{k}{row})/length(cpIdx{row});
    end
end

end