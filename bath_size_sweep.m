%% Parameters
w = 1; 
ww = 5;
gamma = 1;
beta_L = 0.5;
beta_R = 0.5;
mu_L = -0.1;
mu_R = 0.1;
% beta_L = 0.01;
% beta_R = 0.01;
% mu_L = -0.01;
% mu_R = 0.01;
n = 200;
site = n/2;

tic;
%% Evaluation
current_values = zeros(1, n/2-1);
ratio_values = zeros(1, n/2-1);

for nL=1:n/2-1
    nR = nL;
    nW = n-nR-nL;
    ratio_values(nL) = nL/nW;
    current = compute_current_multisite(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, site);
    disp(nL);
    disp(current);
    current_values(nL) = real(current);
end

toc;
%% Plot
figure;
plot(ratio_values, current_values, 'b');
xlabel('Ratio nL/nW');
ylabel('Current');