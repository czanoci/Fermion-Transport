function [ current, V, eigenvalues ] = compute_current( A )

[N, ~] = size(A);
n = N/4;
[eigenvectors, eigenvalues] = eig(A);
eigenvalues = diag(eigenvalues);

V = sort_eigenvalues( eigenvectors, eigenvalues);
% IMPROVE THIS
rounding_const = 1000;
rounded_eigenvalues = floor((abs(real(eigenvalues)) + 0.1/rounding_const)*rounding_const)/rounding_const;
if size(unique(rounded_eigenvalues)) == 1
    num_degen_eigenval = size(rounded_eigenvalues, 1);
else
    [num_degen_eigenval, eigenval] = hist(rounded_eigenvalues, unique(rounded_eigenvalues));
    num_degen_eigenval = sort(num_degen_eigenval, 'descend');
end

% number of different eigenvalues (up to sign)
num_blocks = size(num_degen_eigenval, 2);
processed_eigenvectors = 0;
for block=1:num_blocks
   num_eigenval = num_degen_eigenval(block);
   if mod(num_eigenval, 2) == 1
       disp(eigenvalues);
       disp(rounded_eigenvalues);
   end
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

% D = diag(eigenvalues([6, 5, 1, 4, 8, 2, 3, 7]));
% disp(A*(V.') - (V.')*D);

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

current = -(w1w3+w2w4)/(2*1i);

% w4w3 = 0;
% for m=1:2*n
%    w4w3 = w4w3 + V(2*m, 7)*V(2*m-1, 5) - V(2*m, 8)*V(2*m-1, 6) - 1i*V(2*m, 8)*V(2*m-1, 5) - 1i*V(2*m, 7)*V(2*m-1, 6); 
% end
% 
% w4w3 = w4w3/2;
%disp((1i*w4w3+1)/2);

end

