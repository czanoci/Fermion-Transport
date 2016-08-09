%% Parameters
w = 1; 
ww = 1;
gamma = 1;
mu_L = -0.01;
mu_R = 0.01;
beta_L = 0.1;
beta_R = 0.2;
nL = 4;
nR = 4;
N = 50 - nL - nR;
num_trials = 200;

%% Evaluation
avg_delta_values = zeros(1, N-2);
n_values = zeros(1, N-2);

for i=3:N
    disp(i)
    n_values(i-2) = i + nL + nR;
    delta_values = zeros(1, num_trials);
    for j=1:num_trials
        delta_values(j) = compute_relaxation_rate(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, i, nR, 1.5);
    end
    avg_delta_values(i-2) = mean(delta_values);
end

%% Plot
figure;
scatter(n_values, log(avg_delta_values), 'r');
xlabel('n');
ylabel('log(delta)');
hold on;
X =[ones(length(n_values), 1), n_values.']; 
b = X \ log(avg_delta_values).';
plot(n_values, X*b, 'b');
disp(b);