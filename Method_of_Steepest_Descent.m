function [ x, niters ] = Method_of_Steepest_Descent( A, b, x0 )
% Solves Ax = b using the method of steepest descent
% Inputs:
%   A: matrix A
%   b: vector b
% Outputs:
%   x: solution vector
%   niters: number of iterations performed

% Initialize
x = x0;
r = b - A*x;
niters = 0;

% Iterate until convergence or maximum iterations
while norm(r) > 1e-15 * norm(b)
    p = r;
    q = A * p;
    alpha = dot(p,r) / dot(p,q);
    x = x + alpha*p;
    r = r - alpha*q;
    niters = niters + 1;
end
end