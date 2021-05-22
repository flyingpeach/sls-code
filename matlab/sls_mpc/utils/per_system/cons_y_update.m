function Ys = cons_y_update(Xs, Ys, Zs, rows, row2cp, cpIdx)                

for row = rows
    rIdx = row2cp(row); kIdx = row2cp(cpIdx{row});
    Ys{rIdx} = Ys{rIdx} + Xs{rIdx} - [Zs{kIdx}]';
end

end
             