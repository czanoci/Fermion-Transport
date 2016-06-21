function [ current ] = compute_current_multisite(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR)
%% Write the Hamiltonian of the system in the w_j basis
% See Prosen Eq. (2) 
H = compute_H( w, ww, nL, nW, nR );

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
[N, ~] = size(A);
n = N/4;
[eigenvectors, eigenvalues] = eig(A);
eigenvalues = diag(eigenvalues);

%% Sort the eigenvectors according to their eigenvalues
% The order is beta1, -beta1, beta2, -beta2... where Re(beta1)>=Re(beta2)...
[V, num_degen_eigenval] = sort_eigenvalues_multisite( eigenvectors, eigenvalues);

% number of different eigenvalues (up to sign)
num_blocks = size(num_degen_eigenval, 2);
processed_eigenvectors = 0;
for block=1:num_blocks
   num_eigenval = num_degen_eigenval(block);
   T = zeros(num_eigenval/2, num_eigenval/2);
   for i=1:num_eigenval/2
       for j=1:num_eigenval/2
           T(i, j) = V(processed_eigenvectors+2*i, :)*V(processed_eigenvectors+2*j-1, :).';
       end
   end
   U = (inv(T)).';
   V_copy = V;
   for i=1:num_eigenval/2
       V(processed_eigenvectors+2*i-1, :) = 0;
       for j=1:num_eigenval/2
           V(processed_eigenvectors+2*i-1, :) = V(processed_eigenvectors+2*i-1, :) + U(i, j)*V_copy(processed_eigenvectors+2*j-1, :);
       end
   end
   processed_eigenvectors = processed_eigenvectors + num_eigenval;
end

% Make sure matrix V is normalized correctly
diff = V*(V.') - compute_J(n);
if sum(abs(diff)) > 1E-5;
    disp('Matrix V not normalized correctly');
    disp(sum(abs(diff)));
end

current = -(quadratic_observable(V, 2*nL+nW-1, 2*nL+nW+1)+quadratic_observable(V, 2*nL+nW, 2*nL+nW+2))/(4*1i);
