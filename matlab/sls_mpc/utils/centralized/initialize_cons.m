function [Ys, Zs] = initialize_cons(rCp, cpIdx, nValsCp)                

Ys = cell(nValsCp, 1);
Zs = cell(nValsCp, 1);

for i=1:length(rCp)
    for row = rCp{i}
        Zs{row} = 0;
        for k = cpIdx{row}
            Ys{row}{k} = 0;
        end
    end
end

end
             