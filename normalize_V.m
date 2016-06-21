function [ norm_V ] = normalize_V( V, num_degen_eigenval )
% Take a linear combination of eigenvectors so that the entries in V are normalized as described in Eq. (30)

norm_V = zeros(size(V));

% Number of blocks in which the absolute value of real part of the
% eigenvalue is the same. 
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
   for i=1:num_eigenval/2
       norm_V(processed_eigenvectors+2*i-1, :) = 0;
       norm_V(processed_eigenvectors+2*i, :) = V(processed_eigenvectors+2*i, :);
       for j=1:num_eigenval/2
           norm_V(processed_eigenvectors+2*i-1, :) = norm_V(processed_eigenvectors+2*i-1, :) + U(i, j)*V(processed_eigenvectors+2*j-1, :);
       end
   end
   processed_eigenvectors = processed_eigenvectors + num_eigenval;
end

end

