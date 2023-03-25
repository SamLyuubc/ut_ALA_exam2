function [ x, niters ] = CG( A, b, x0 )
% Initialize
x = x0;
r = b;
niters = 0;
p_old = r;

% Iterate until convergence or maximum iterations
while norm(r) > 1e-15 * norm(b)
    if niters == 0
        p = r;
    else
        ga = dot(r,r)/dot(r_old, r_old);
        p_old = p;
        p = r + ga * p_old;
    end
    alpha = dot(p,r) / dot(p,A*p);
    x = x + alpha*p;
    r_old = r;
    r = r - alpha*A*p;
    niters = niters + 1;
end
end