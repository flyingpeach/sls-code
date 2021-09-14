function [A_tilde, b_tilde] = project_set(A, b, len_x)
% Find H, h such that the set
%   x s.t. A_tilde * x <= b_tilde
% is equivalent to
%   x s.t. \exists w s.t. A*[w; x] <= b
% where x is vector of length len_x

len_wx = size(A, 2);
len_w  = len_wx - len_x;

for i=1:len_w
    [A_tilde, b_tilde] = fourier_motzkin(A, b);
    
    % Reduce redundant constraints via mpt3 toolbox
    P = Polyhedron(A_tilde, b_tilde);
    P.minHRep();

    A_tilde = P.A; b_tilde = P.b;
    A = A_tilde; b = b_tilde;
end

end