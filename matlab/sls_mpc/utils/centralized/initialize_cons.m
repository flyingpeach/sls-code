function [Ys, Zs] = initialize_cons(rCp, cpIdx, nCp)                

Ys = cell(nCp, 1);
Zs = cell(nCp, 1);

for i=1:length(rCp)
    for row = rCp{i}
        Zs{row} = 0;
        for k = cpIdx{row}
            Ys{row}{k} = 0;
        end
    end
end

end
             