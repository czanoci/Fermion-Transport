function [ corr_matrix ] = compute_corr_matrix( V )
% Compute the correlation matrix M with entries:
% M_jk = < c_j^\dagger c_k > where j, k belong to B
% Correlation matrix is written in the c_j basis
% B is a subset of the chain, represented by a vector of indices
% if the index i is on the list, then site i belongs to B

% n = size(B, 2);
% corr_matrix = zeros(n, n); % mp
% for j = 1:n
%     for k = 1:n
%         % use index to get the correct site
%         x = B(j);
%         y = B(k);
%         corr_matrix(j, k) = two_point_corr( x, y, V );
%     end
% end

N = size(V, 2);
n = N/4;
corr_matrix = zeros(n, n);
for i = 1:n
    for j = 1:n
        corr_matrix(i, j) = two_point_corr( i, j, V );
    end
end
end

