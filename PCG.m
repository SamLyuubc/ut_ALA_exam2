function [ x, niters ] = PCG( A, b, x0 )
% Initialize
x = x0;
r = b;
niters = 0;

L = sparse( ichol(sparse(A), struct('type','ict','droptol',1e-4)));
M = L * L';
z = r;

% Iterate until convergence or maximum iterations
while norm(r) > 1e-15 * norm(b)
    z_old = z;
    z = M\r;
    if niters == 0
        p = z;
    else
        ga = dot(r, z)/ dot(r_old, z_old);
        p_old = p;
        p = z + ga * p_old;
    end
    q = A * p;
    alpha = dot(p,r) / dot(p,q);
    x = x + alpha*p;
    r_old = r;
    r = r - alpha*q;
    niters = niters + 1;
end