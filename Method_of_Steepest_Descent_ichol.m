function [ x, niters ] = Method_of_Steepest_Descent_ichol( A, b, x0 )
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
L = sparse( ichol(sparse(A), struct('type','ict','droptol',1e-4)));
M = L * L';

% Iterate until convergence or maximum iterations
while norm(r) > 1e-10
    p = M\r;
    q = A * p;
    alpha = dot(p,r) / dot(p,q);
    x = x + alpha*p;
    r = r - alpha*q;
    niters = niters + 1;
end
end