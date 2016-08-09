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

for i=1:4*n
    aug_eigenvectors(i, 1) = real(eigenvalues(i));
    aug_eigenvectors(i, 2) = imag(eigenvalues(i));
    aug_eigenvectors(i, 3:end) = eigenvectors(:, i);
end

dot_products = (eigenvectors.')*eigenvectors;
threshold = 0.1;
dot_products( abs(dot_products) >= threshold ) = 1;
dot_products( abs(dot_products) < threshold ) = 0;
seen = [];
groups = {};
num_groups = 0;

for i = 1:size(dot_products, 1)
   if ismember(i, seen) == 0
       num_groups = num_groups + 1;
       row = find( dot_products(i, :) == 1 );
       col = find( dot_products(:, row(1)) == 1).';
       seen = cat(2, seen, row);
       seen = cat(2, seen, col);
       group = aug_eigenvectors([row, col], :);
       groups{num_groups} = group;
   end
end

unique_real_parts = zeros(1, num_groups);

for j = 1:num_groups
   E = groups{j};
   unique_real_parts(j) = min(E(:, 1));
   pos_real_E = E(E(:, 1) >= 0, :);
   neg_real_E = E(E(:, 1) < 0, :);
   pos_pos_E = pos_real_E(pos_real_E(:, 2) >= 0, 3:end);
   pos_neg_E = pos_real_E(pos_real_E(:, 2) < 0, 3:end);
   neg_pos_E = neg_real_E(neg_real_E(:, 2) >= 0, 3:end);
   neg_neg_E = neg_real_E(neg_real_E(:, 2) < 0, 3:end);
   ordered_E = [pos_pos_E; neg_neg_E; pos_neg_E; neg_pos_E];
   groups{j} = ordered_E;
end

[~, index] = sort(unique_real_parts, 'descend');
V = [];
num_degen_eigenval = zeros(1, num_groups);

for j = 1:num_groups
    E = groups{index(j)};
    num_degen_eigenval(j) = size(E, 1);
    V = cat(1, V, E);
end

end

