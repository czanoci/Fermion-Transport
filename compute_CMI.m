function [ CMI ] = compute_CMI( B, C, D, corr_matrix, S_BCD )
% Compute conditional mutual information between B, C, and D
% CMI(B, C, D) = S(BC) + S(CD) - S(C) - S(BCD)

corr_matrix_C = corr_matrix(C, C);%compute_corr_matrix( C, V );
corr_matrix_BC = corr_matrix(union(B, C), union(B, C));%compute_corr_matrix( union(B, C), V );
corr_matrix_CD = corr_matrix(union(C, D), union(C, D));%compute_corr_matrix( union(C, D), V );
%corr_matrix_BCD = compute_corr_matrix( union(union(B, C), D), V );
S_C = entropy(corr_matrix_C);
S_BC = entropy(corr_matrix_BC);
S_CD = entropy(corr_matrix_CD);
%S_BCD = entropy(corr_matrix_BCD);
CMI = S_BC + S_CD - S_C - S_BCD;

if real(CMI) < 0
    disp('WARNING: Conditional mutual information is negative');
    disp(real(CMI));
end

end

