w = 10;
ww = 10;
nL = 20;
nW = 0;
nR = 0;
gamma = 2;
beta_L = 0.1;
beta_R = 0;
mu_L = 0.5;
mu_R = 0;

H = compute_H( w, ww, nL, nW, nR );
L = compute_L( w, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR );
M = compute_M(L);
A = compute_A(H, M);

if A ~= -A.'
    disp('Matrix A not anti-symmetric');
end

[N, ~] = size(A);
n = N/4;
[eigenvectors, eigenvalues] = eig(A);
eigenvalues = diag(eigenvalues);

[V, num_degen_eigenval] = sort_eigenvalues_old( eigenvectors, eigenvalues);
if sum(mod(num_degen_eigenval, 2)) ~= 0
    disp('num_degen_eigenval has odd entries');
end

V = normalize_V( V, num_degen_eigenval );

% Make sure matrix V is normalized correctly
diff = V*(V.') - compute_J(n);
if sum(abs(diff(:))) > 1E-3;
    disp('Matrix V not normalized correctly');
    disp(max(abs(diff(:))));
    disp(sum(abs(diff(:))));
end

[occupation_nums, energies] = occupation_num(V, w, 0, nL);
th_occupation_nums = 1./(exp(beta_L*(energies-mu_L))+1);

figure
scatter(energies, th_occupation_nums, 'filled', 'r');
hold on;
scatter(energies, occupation_nums, 'b');
xlabel('Energy');
ylabel('Occupation Number');
title(['Occupation number vs energy for beta = ', num2str(beta_L)]);
legend('theoretically predicted', 'numerically computed');
