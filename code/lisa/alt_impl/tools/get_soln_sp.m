function nullSp = get_soln_sp(F2, tol)
% Calculate size of solution space for each Tc
% Outputs
%     nullSp        : the nullspace of F2
% Inputs
%     F1s           : LTISystem containing system matrices
%     slsParams     : SLSParams containing parameters
%     slsOuts       : contains info from SLS (original R, M)
%     Tc            : length of the approximate solution
%     tol           : tolerance to use for rank calculation

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




    

