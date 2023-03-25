function [ A ] = Create_Poisson_problem_A( N )

  % Create the archtypical matrix A for an N x N Poisson problem (5-point
  % stencil.

  % Set the diagonal

  
  % Set the entries of the first sub and super diagonals

  
  % Set the other off-diagonal entries

  A = sparse(N*N, N*N);   % initialize the matrix as a sparse matrix

  % Set the diagonal entries
  for i = 1:N
      for j = 1:N
        k = (i-1)*N + j;
        A(k,k) = 4;
      end
  end

  % Set the entries of the first sub and super diagonals
  for i = 1:N-1
      for j = 1:N-1
          k = (i-1)*N + j;
          A(k,k+1) = -1;
          A(k+1,k) = -1;
      end
  end

  % Set the other off-diagonal entries
  for i = 1:N-2
      for j = 1:N-2
          k = (i-1)*N + j;
          A(k,k+N) = -1;
          A(k+N,k) = -1;
      end
  end
end