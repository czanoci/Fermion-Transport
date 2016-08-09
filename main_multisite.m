%% Parameters
w = 1; 
ww = 1;
gamma = 1;
nL = 5;
nW = 20;
nR = 5;
site = floor(nL+nW/2);
mu_R = 0.1;
mu_L = 0.1;

tic;
%% Evaluation
delta_beta_values = zeros(1, 5);
current_values = zeros(1, 5);
kappa_values = zeros(1, 10);
beta_values = zeros(1, 10);

for i=1:10
    beta = 0.01*i;
    beta_values(i) = beta;
    for j=1:5
        beta_L = (1-0.05*j)*beta;
        beta_R = (1+0.05*j)*beta;
        delta_beta_values(j) = beta_R - beta_L;
        current_values(j) = real(compute_current_multisite(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, site, 0));
    end
    slope = current_values / delta_beta_values;
    kappa_values(i) = slope;
    disp(i);
end

% w = 1, ww =1, nL=nR=5, nW=10, beta = 0.01
%mathematica_sigma_values = [0.00122016, 0.00244005, 0.00365941, 0.00487796, 0.00609546, 0.00731162, 0.00852619, 0.0097389, 0.0109495, 0.0121577];

% w = 1, ww =1, nL=nR=6, nW=14, beta = 0.01
%mathematica_sigma_values = [0.000991833, 0.00198347, 0.00297473, 0.0039654, 0.0049553, 0.00594424, 0.00693202, 0.00791846, 0.00890335, 0.00988652];

% w = 15, ww = 5, nL=nR=6, nW = 14, beta=0.1
%mathematica_sigma_values = [ 0.00132812, 0.001555, 0.00167423, 0.001418, 0.000599181, -0.000612653, 1.61188*1E-7, 1.846*1E-7, 2.08167*1E-7, 2.31908*1E-7];

% w = 15, ww = 5, nL=nR=6, nW = 14, beta=0.01
%mathematica_sigma_values = [ 0.000184101, 0.000363236, 0.000533031, 0.000690132, 0.00083244, 0.000959171, 0.00107075, 0.00116843, 0.00125379, 0.00132812];
toc;
%% Plot
figure;
%plot(beta_values, mathematica_sigma_values, 'b');
xlabel('Beta');
ylabel('Conductivity kappa');
hold on;
plot(beta_values, kappa_values, 'r');
%legend('Mathematica computation', 'Prosen method');