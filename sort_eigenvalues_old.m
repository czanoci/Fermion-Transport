function [ V, num_degen_eigenval ] = sort_eigenvalues( eigenvectors, eigenvalues)
% Sorts eigenvectors such that their eigenvalues are arranged in the form
% -a-b*i = -beta1
% -a+b*i = -beta2
% a-b*i = beta2
% a+b*i = beta1
% then the jth and (2n-j)th eigenvectors are paired (put next to each other
% in V)

% eigenvectors augmented with a and b at the beginning
% a is the real part of eigenvalue, b is the imaginary part
aug_eigenvectors = zeros(size(eigenvectors, 2), size(eigenvectors, 1) + 2);
n = size(eigenvalues, 1) / 4;

% Some eigenvalues should have the same real part, but due to numerical
% precision, they may be slightly different, e.g. 1.83249999 and 1.8325000
% In order to sort the eigenvalues correctly, I only use 3 significant
% digits when comparing them. This is done using digits and vpa functions
digits(3);
for i=1:4*n
    aug_eigenvectors(i, 1) = vpa(real(eigenvalues(i)));
    aug_eigenvectors(i, 2) = vpa(imag(eigenvalues(i)));
    aug_eigenvectors(i, 3:end) = eigenvectors(:, i);
end

% Now sort the eigenvectors based on the real part of the eigenvalues (col 1).
% Break ties based on the imaginary part of the eigenvalues. (col 2)
[~, sortedIdx] = sortrows(real(aug_eigenvectors), [1, 2]);
aug_eigenvectors = aug_eigenvectors(sortedIdx, :);

% Compute the number of eigenvalues that have the same real part (up to sign)
%disp(abs(aug_eigenvectors(:, 1)));
if size(unique(abs(aug_eigenvectors(:, 1))), 1) > 1
    [num_degen_eigenval, ~] = hist(abs(aug_eigenvectors(:, 1)), unique(abs(aug_eigenvectors(:, 1))));
    num_degen_eigenval = fliplr(num_degen_eigenval);
else 
    num_degen_eigenval = size(aug_eigenvectors, 1);
end

%disp(num_degen_eigenval);
V = zeros(4*n, 4*n);
% Pair the eigenvectors in the entries of V
for i=1:2*n
   V(2*i-1, :) = aug_eigenvectors(4*n+1-i, 3:end); % pos eig
   %disp(aug_eigenvectors(4*n+1-i, 1) + 1i*aug_eigenvectors(4*n+1-i, 2));
   V(2*i, :) = aug_eigenvectors(i, 3:end); % neg eig
   %disp(aug_eigenvectors(i, 1) + 1i*aug_eigenvectors(i, 2));
end

end

