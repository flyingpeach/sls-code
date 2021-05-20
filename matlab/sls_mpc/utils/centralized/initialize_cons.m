function [Ys, Zs] = initialize_cons(rCp, cpIdx, row2cp, nCp)                

Ys = cell(nCp, 1);
Zs = cell(nCp, 1);

for i=1:length(rCp)
    for row = rCp{i}
        rIdx = row2cp(row); kIdx = row2cp(cpIdx{row});
        Zs{rIdx} = 0;
        for k = kIdx
            Ys{rIdx}{k} = 0;
        end
    end
end

end
             