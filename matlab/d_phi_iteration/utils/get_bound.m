function beta = get_bound(D, M, stabType)

if stabType == StabType.L1
    beta = max(sum(D*M/D, 2));
elseif stabType == StabType.LInf
    beta = max(sum(D*M/D, 1));
else % default to nu
    beta = max(max(D*M/D));
end

end