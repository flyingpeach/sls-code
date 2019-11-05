function nullSp = get_nullsp(F2, tol)
% Calculates the nullspace of F2 with the specified tolerance
% i.e. right singular vectors corresponding to singular values lower than
% tol

[~, S, V] = svd(F2); % finding nullspace with tolerance

idxZer = 0;
for i=1:size(F2, 2)
    if S(i,i) <= tol
        idxZer = i;
        break;
    end
end

if idxZer == 0
    nullSp = [];
else
    nullSp = V(:, idxZer:size(F2, 2));
end