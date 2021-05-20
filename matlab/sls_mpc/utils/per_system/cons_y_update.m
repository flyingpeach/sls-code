function Ys = cons_y_update(Xs, Ys, Zs, rows, row2cp, cpIdx)                

for row = rows
    rIdx = row2cp(row); kIdx = row2cp(cpIdx{row});
    for k = kIdx % neighbors
        Ys{rIdx}{k} = Ys{rIdx}{k} + Xs{rIdx}(k) - Zs{k};
    end
end

end
             