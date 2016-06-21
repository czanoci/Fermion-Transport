%% Parameters
w = 1; 
ww = 1;
gamma = 1;
nL = 5;
nW = 20;
nR = 5;

%% Evaluation
mu_values = zeros(1, 5);
current_values = zeros(1, 5);
sigma_values = zeros(1, 10);
beta_values = zeros(1, 10);

for i=1:10
    beta_L = 0.01*i;
    beta_R = 0.01*i;
    beta_values(i) = 0.01*i;
    
    for j=1:5
        mu_L = -0.01*j;
        mu_R = 0.01*j;
        mu_values(j) = mu_R-mu_L;
        current_values(j) = -2*real(compute_current_multisite(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR));
    end
    slope = current_values / mu_values;
    sigma_values(i) = slope;
end

% ww =1, nL=nR=5, nW=10
%mathematica_sigma_values = [0.00122016, 0.00244005, 0.00365941, 0.00487796, 0.00609546, 0.00731162, 0.00852619, 0.0097389, 0.0109495, 0.0121577];

% ww =1, nL=nR=6, nW=14
mathematica_sigma_values = [0.000991833, 0.00198347, 0.00297473, 0.0039654, 0.0049553, 0.00594424, 0.00693202, 0.00791846, 0.00890335, 0.00988652];
%% Plot
figure;
% plot(beta_values, mathematica_sigma_values, 'b');
xlabel('Beta');
ylabel('Conductivity sigma');
% hold on;
plot(beta_values, sigma_values, 'r');
% legend('Mathematica computation', 'Prosen method');