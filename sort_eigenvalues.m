function [ V ] = sort_eigenvalues( eigenvectors, eigenvalues)
% Sorts eigenvectors such that their eigenvalues are arranged in the form
% -a-b*i
% -a+b*i
% a-b*i
% a+b*i
% then the jth and (2n-j)th eigenvectors are paired (put next to each other
% in V)

% eigenvectors augmented with a and b at the beginning
aug_eigenvectors = zeros(size(eigenvectors, 2), size(eigenvectors, 1) + 2);
n = size(eigenvalues, 1) / 4;

for i=1:4*n
    aug_eigenvectors(i, 1) = real(eigenvalues(i));
    aug_eigenvectors(i, 2) = imag(eigenvalues(i));
    aug_eigenvectors(i, 3:end) = eigenvectors(:, i);
end
[~, sortedIdx] = sortrows(real(aug_eigenvectors), [1, 2]);
aug_eigenvectors = aug_eigenvectors(sortedIdx, :);

V = zeros(4*n, 4*n);

for i=1:2*n
   V(2*i-1, :) = aug_eigenvectors(4*n+1-i, 3:end); % pos eig
   V(2*i, :) = aug_eigenvectors(i, 3:end); % neg eig
end

% V(1, :) = aug_eigenvectors(7, 3:end);
% V(2, :) = aug_eigenvectors(1, 3:end);
% V(3, :) = aug_eigenvectors(8, 3:end);
% V(4, :) = aug_eigenvectors(2, 3:end);
% V(5, :) = aug_eigenvectors(6, 3:end);
% V(6, :) = aug_eigenvectors(4, 3:end);
% V(7, :) = aug_eigenvectors(5, 3:end);
% V(8, :) = aug_eigenvectors(3, 3:end);

end

