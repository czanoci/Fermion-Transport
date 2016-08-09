%% Parameters
w = 1; 
ww = 3;
gamma = 1;
beta_L = 0.01;
beta_R = 0.01;
mu_L = -0.01;
mu_R = 0.01;
nL = 20;
nR = 60;
nW = 20;
n = nL+nW+nR;
subset_size = 10;
B = [(n - subset_size)/2 : (n + subset_size)/2 - 1];

tic;
%% Evaluation
MI_values = zeros(1, n - subset_size + 1);
position_values = zeros(1, n - subset_size + 1);

for pos=1 : n - subset_size + 1
    position_values(pos) = pos;
    C = [pos : pos + subset_size - 1];
    MI = mutual_info(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, B, C);
    disp(MI);
    MI_values(pos) = real(MI);
end

toc;
%% Plot
figure;
plot(position_values, MI_values, 'b');
xlabel('Measurement site');
ylabel('Current');