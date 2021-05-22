function Ys = cons_y_update(Xs, Ys, Zs, rows, row2cp, cpIdx)                

for row = rows
    rIdx = row2cp(row); kIdx = row2cp(cpIdx{row});
    for j = 1:length(kIdx) % neighbors        
        Ys{rIdx}{j} = Ys{rIdx}{j} + Xs{rIdx}(j) - Zs{kIdx(j)};
    end
end

end
             