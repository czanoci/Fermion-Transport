%% Parameters
w = 1; 
ww = 0.1;
gamma_values = [0.01, 0.1, 1];
beta = 5;
beta_L = beta;
beta_R = beta;
mu_L = 0;
mu_R = 0;
nL = 160;
nR = 160;
nW = 160;
n = nL+nW+nR;
disorder = 0;


tic;
%% Evaluation
figure;
for i = 1:size(gamma_values, 2)
    gamma = gamma_values(i);
    [V, error] = mutual_info(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, disorder);
    [occupation_nums, energies] = occupation_num(V, w, nL, nW);
    th_occupation_nums = 1./(exp(beta*(energies))+1);
    scatter(energies, occupation_nums);
    hold on;
end

toc;
%% Plot

scatter(energies, th_occupation_nums, 'filled', 'y');
xlabel('Energy');
ylabel('Occupation Number');
% title(['Occupation number vs energy for beta = ', num2str(beta_L), ', ww = ', num2str(ww), ', gamma =', num2str(gamma)]);
% legend('Fermi distribution', 'Computed numerically');

% figure;
% scatter(n_values, occupation_values, 'b');
