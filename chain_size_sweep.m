%% Parameters
w = 1; 
ww = 2;
gamma = 1;
beta_L = 0.5;
beta_R = 0.5;
mu_L = -0.1;
mu_R = 0.1;
n_max = 400;
ratio_nL_to_n = 1/4;

tic;
%% Evaluation
current_values = zeros(1, 100);
n_values = zeros(1, 100);

for nL=1:n_max*ratio_nL_to_n
    nR = nL;
    n = nL/ratio_nL_to_n;
    n_values(nL) = n;
    nW = n - nL - nR;
    site = n/2;
    current = compute_current_multisite(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, site);
    disp(nL);
    disp(current);
    current_values(nL) = real(current);
    n_values(nL) = n;
end

toc;
%% Plot
figure;
plot(n_values, current_values, 'b');
xlabel('Total chain size n');
ylabel('Current');