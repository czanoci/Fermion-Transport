function [ S ] = entropy( corr_matrix )
% Compute the Shannon entropy from correlation matrix

I = eye(size(corr_matrix));
S = -trace(corr_matrix * logm(corr_matrix) + (I - corr_matrix) * logm(I - corr_matrix));

end

