function Zs = cons_z_update(Xs, Ys, Zs, rows, row2cp, cpIdx)                

for row = rows
    rIdx = row2cp(row);
    Zs{rIdx} = 0;
    for k = cpIdx{row} % neighbors
        j = find(cpIdx{k} == row); % the index where your neighbor stores YOUR value
        kIdx = row2cp(k);
        Zs{rIdx} = Zs{rIdx} + (Xs{kIdx}(j)+Ys{kIdx}(j))/length(cpIdx{row});
    end
end

end