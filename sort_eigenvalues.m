function [ V ] = sort_eigenvalues( eigenvectors, eigenvalues )
% Sorts eigenvectors such that their eigenvalues are arranged in the form
% -a-b*i
% -a+b*i
% a-b*i
% a+b*i
% then the jth and (2n-j)th eigenvectors are paired (put next to each other
% in V)

% eigenvectors augmented with a and b at the beginning
aug_eigenvectors = zeros(size(eigenvectors, 2), size(eigenvectors, 1) + 2);
n = size(eigenvalues, 1) / 2;

for i=1:2*n
    aug_eigenvectors(i, 1) = real(eigenvalues(i));
    aug_eigenvectors(i, 2) = imag(eigenvalues(i));
    aug_eigenvectors(i, 3:end) = eigenvectors(:, i);
end

aug_eigenvectors = sortrows(aug_eigenvectors, [1, 2]);
%disp(aug_eigenvectors);
V = zeros(2*n, 2*n);

for i=1:n
   V(2*i-1, :) = aug_eigenvectors(i, 3:end); % neg eig
   V(2*i, :) = aug_eigenvectors(2*n+1-i, 3:end); % pos eig
end
end

