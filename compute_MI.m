function [ MI ] = compute_MI( B, C, corr_matrix )
% Compute mutual information between B and C
% MI(B, C) = S(C) + S(B) - S(BC)

corr_matrix_B = corr_matrix(B, B);%compute_corr_matrix( B, V );
corr_matrix_C = corr_matrix(C, C);%compute_corr_matrix( C, V );
corr_matrix_BC = corr_matrix(union(B, C), union(B, C));%compute_corr_matrix( union(B, C), V );
S_B = entropy(corr_matrix_B);
S_C = entropy(corr_matrix_C);
S_BC = entropy(corr_matrix_BC);
MI = S_B + S_C - S_BC;

if real(MI) < 0
    disp('WARNING: Mutual information is negative');
    disp(real(MI));
end

end

