function [ delta ] = compute_relaxation_rate(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, disorder)
%% Write the Hamiltonian of the system in the w_j basis
% See Prosen Eq. (2) 
H = compute_H( w, ww, nL, nW, nR )+ add_disorder( disorder, nL, nW, nR );

%% Matrix of Lindblad operators
% where each row is a bath operator L_mu written in the w_j basis 
% L_j's are multiplied 1/sqrt(2) to account for diff. in Lindblad eq. in
% Prosen and transport papers. The eq. in Prosen has an extra factor of 2!
% See Prosen Eq. (3)
L = compute_L( w, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR );

%% Compute the matrix M defined in Prosen Eq. (23)
M = compute_M(L);

%% Compute the matrix A defined in Prosen Eq. (27)
A = compute_A(H, M);

% Make sure matrix A is anti-symmetric
if A ~= -A.'
    disp('Matrix A not anti-symmetric');
end

%% Compute eigenvalues and eigenvectors
[~, eigenvalues] = eig(A);
eigenvalues = diag(eigenvalues);

%% Compute rate of relaxation
real_eigenvalues = real(eigenvalues);
pos_real_eigenvalues = real_eigenvalues(real_eigenvalues >= 0);
delta = 2*min(pos_real_eigenvalues);

end