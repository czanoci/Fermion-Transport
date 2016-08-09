%% Parameters
w = 1; 
ww = 1;
gamma = 1;
beta_L = 2;
beta_R = 3;
mu_L = -0.01;
mu_R = 0.01;
nL = 5;
nR = 5;
nW = 40;
site = 25;
V_values = [0, 0.01*w, 0.02*w, 0.05*w, 0.1*w, 0.2*w, 0.5*w, w, 1.5*w, 2*w];
num_trials = 200;

tic;
%% Evaluation
avg_current_values = zeros(1, 10);

for i=1:10
    disp(i)
    current_values = zeros(1, num_trials);
    for j=1:num_trials
        current = compute_current_multisite(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, site, V_values(i));
        current_values (j) = real(current);
    end
    avg_current_values(i) = mean(current_values)
end
disp(avg_current_values);

toc;
%% Plot
figure;
plot(V_values, avg_current_values, 'b');
xlabel('Disorder max strength, V0');
ylabel('Average current in the middle of chain');