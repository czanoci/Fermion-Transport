function [ current ] = compute_current( A )
% Given the matrix A, compute the current using Theorem 3 from Prosen's paper

%% Compute eigenvalues and eigenvectors
[N, ~] = size(A);
n = N/4;
[eigenvectors, eigenvalues] = eig(A);
eigenvalues = diag(eigenvalues);

%% Sort the eigenvectors according to their eigenvalues
% The order is beta1, -beta1, beta2, -beta2... where Re(beta1)>=Re(beta2)...
V = sort_eigenvalues( eigenvectors, eigenvalues);

%% Take a linear combination of eigenvectors so that the entries in V are normalized as described in Eq. (30)

% Number of blocks in which the absolute value of real part of the
% eigenvalue is the same. In a 2-site model, there are 2 blocks with 4
% eigenvectors each
num_degen_eigenval = [4, 4];

% number of different eigenvalues (up to sign)
num_blocks = size(num_degen_eigenval, 2);
processed_eigenvectors = 0;
for block=1:num_blocks
   num_eigenval = num_degen_eigenval(block);
   T = zeros(num_eigenval/2, num_eigenval/2);
   for i=1:num_eigenval/2
       for j=1:num_eigenval/2
           T(i, j) = V(processed_eigenvectors+2*i, :)*V(processed_eigenvectors+2*j-1, :).';
       end
   end
   U = (inv(T)).';
   V_copy = V;
   for i=1:num_eigenval/2
       V(processed_eigenvectors+2*i-1, :) = 0;
       for j=1:num_eigenval/2
           V(processed_eigenvectors+2*i-1, :) = V(processed_eigenvectors+2*i-1, :) + U(i, j)*V_copy(processed_eigenvectors+2*j-1, :);
       end
   end
   processed_eigenvectors = processed_eigenvectors + num_eigenval;
end

%% Use Theorem 3 / Eq. (47) to compute the quadratic observables wiwj
w1w3 = 0;
for m=1:2*n
   w1w3 = w1w3 + V(2*m, 1)*V(2*m-1, 5) - V(2*m, 2)*V(2*m-1, 6) - 1i*V(2*m, 2)*V(2*m-1, 5) - 1i*V(2*m, 1)*V(2*m-1, 6); 
end

w1w3 = w1w3/2;

w2w4 = 0;
for m=1:2*n
   w2w4 = w2w4 + V(2*m, 3)*V(2*m-1, 7) - V(2*m, 4)*V(2*m-1, 8) - 1i*V(2*m, 4)*V(2*m-1, 7) - 1i*V(2*m, 3)*V(2*m-1, 8); 
end

w2w4 = w2w4/2;

%% Finally, compute the current
% there is an extra factor of 1/2 to match the output of current_method1
current = -(w1w3+w2w4)/(4*1i);

end

