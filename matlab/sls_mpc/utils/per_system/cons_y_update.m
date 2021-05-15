function Ys = cons_y_update(Xs, Ys, Zs, rows, cpIdx)                

for row = rows
    for k = cpIdx{row} % neighbors
        Ys{row}{k} = Ys{row}{k} + Xs{row}(k) - Zs{k};
    end
end

end
             