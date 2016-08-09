function [ occupation_nums, energies ] = occupation_num( V, w, nL, nW )
% compute occupation number using single particle energies/eigenfunctions
% of the middle part of the wire

occupation_nums = zeros(1, nW);

hW = bath_hamiltonian( w, nW );

[eigenvectors_W, eigenvalues_W] = eig(hW);
energies = diag(eigenvalues_W).';
% eigenvectors are columns of eigenvectors_W

for i = nL+1:nL+nW
    for j = nL+1:nL+nW
        ci_dagger_cj = two_point_corr( i, j, V );
        for k = 1:nW
            psi_k = eigenvectors_W(:, k);
            occupation_nums(k) = occupation_nums(k) + psi_k(i-nL)*conj(psi_k(j-nL))*ci_dagger_cj;
        end
    end
end

% for n = 1:nW
%     for x = nL+1:nL+nW
%         for y = nL+1:nL+nW
%             occupation_nums(1, n) = occupation_nums(1, n) + eigenvectors_W(x, n)*conj(eigenvectors_W(y, n))*two_point_corr(x, y, V);
%         end
%     end
% end

end

