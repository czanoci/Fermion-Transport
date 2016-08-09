function [ V, error ] = mutual_info(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, dis)

error = 0;
%% Write the Hamiltonian of the system in the w_j basis
% See Prosen Eq. (2) 
H = compute_H( w, ww, nL, nW, nR ) + add_disorder( dis, nL, nW, nR );% + aubry_andre_model( dis, nL, nW, nR, 0.02*sqrt(2), 0 ); %add_disorder( dis, nL, nW, nR ); 

%% Matrix of Lindblad operators
% where each row is a bath operator L_mu written in the w_j basis 
% L_j's are multiplied 1/sqrt(2) to account for diff. in Lindblad eq. in
% Prosen and transport papers. The eq. in Prosen has an extra factor of 2!
% See Prosen Eq. (3)
L = compute_L( w, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR );

%% Compute the matrix M defined in Prosen Eq. (23)
M = compute_M(L);

%% Compute the matrix A defined in Prosen Eq. (27)
A = compute_A(H, M); % mp

% Make sure matrix A is anti-symmetric
if A ~= -A.'
    disp('Matrix A not anti-symmetric');
end

%% Compute eigenvalues and eigenvectors
[N, ~] = size(A);
n = N/4;
[eigenvectors, eigenvalues] = eig(A);
eigenvalues = diag(eigenvalues);

%% Sort the eigenvectors according to their eigenvalues
% The order is beta1, -beta1, beta2, -beta2... where Re(beta1)>=Re(beta2)...
[V, num_degen_eigenval] = sort_eigenvalues_old( eigenvectors, eigenvalues);
if sum(mod(num_degen_eigenval, 2)) ~= 0
    disp('num_degen_eigenval has odd entries');
    error = 1;
    return;
end

%% Normalize eigenvectors in V so that they satisfy Eq. (30)
V = normalize_V( V, num_degen_eigenval );

% Make sure matrix V is normalized correctly
diff = V*(V.') - compute_J(n);
if sum(abs(diff(:))) > 1E-3;
    disp('Matrix V not normalized correctly');
    disp(max(abs(diff(:))));
    disp(sum(abs(diff(:))));
    error = 1;
end

