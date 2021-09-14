function [A_tilde, b_tilde] = fourier_motzkin(A, b, pos, varargin)
% Find A_tilde, b_tilde such that the set 
%   y_2 s.t. A_tilde * y_2 <= b_tilde
% is equivalent to
%   y_2 s.t. \exists y_1 s.t. A*[y1; y2] <= b
% where y_1 is a scalar

if nargin > 2 % If the user defines the position of y1
	A = [A(:,pos) A(:,1:pos-1) A(:,pos+1:end)];
end

eps = 1e-6; % Effectively zero
n   = size(A, 2);

I_pos  = find(A(:,1) > eps)';
I_neg  = find(A(:,1) < -eps)';
I_zero = find(abs(A(:,1)) < eps)';

A_prime = A;
b_prime = b;
for i = [I_pos I_neg]
    A_prime(i, :) = A(i, :) ./ abs(A(i, 1));
    b_prime(i)    = b(i)    ./ abs(A(i, 1));
end

zeye = [zeros(1, n-1);  eye(n-1)];

A_tilde = [];
b_tilde = [];
if ~isempty(I_zero)
    A_tilde = A(I_zero, :) * zeye;
    b_tilde = b(I_zero);
end

% Add rows corresponding to I_pos and I_neg
for i = I_pos
    for k = I_neg
        A_tilde = [A_tilde; (A_prime(i, :) + A_prime(k, :)) * zeye];
        b_tilde = [b_tilde; b_prime(i) + b_prime(k)];        
    end
end

end
