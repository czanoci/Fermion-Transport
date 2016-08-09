%% Parameters
w = 1; 
ww = 1;
gamma = 1;
% beta_L = 0.1;
% beta_R = 0.1;
% mu_L = -0.01;
% mu_R = 0.01;
beta_values = [1];
mu_values = [0.1];
nL = 20;
nR = 20;
nW = 60;
n = nL+nW+nR;
disorder_values = [0, 0.1];
num_trials = size(disorder_values, 2);

tic;
%% Evaluation
for j=1:size(mu_values, 2)
    for i=1:size(beta_values, 2)
        current_values = zeros(num_trials, 1);
        MI_values = zeros(num_trials, n/2 - 1);
        CMI_values = zeros(num_trials, n/2 - 1);
        beta_L = beta_values(i);
        beta_R = beta_values(i);
        mu_L = -mu_values(j);
        mu_R = mu_values(j);
        for trial=1:num_trials
            disp(trial);
            disorder = disorder_values(trial);
            error = 1;
            while error > 0
                [V, error] = mutual_info(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, disorder);
                disp('Here');
            end
            corr_matrix = compute_corr_matrix( V );
            S_BCD = entropy(corr_matrix);
            for sizeC=1 : n/2-1
                disp(sizeC);
                B = [1 : (n/2 - sizeC)];
                C = [(n/2 - sizeC + 1) : (n/2 + sizeC)];
                D = [(n/2 + sizeC + 1) : n];
                MI = compute_MI( B, D, corr_matrix );
                MI_values(trial, sizeC) = real(MI);
                CMI = compute_CMI( B, C, D, corr_matrix, S_BCD);
                CMI_values(trial, sizeC) = real(CMI);
            end
            current_values(trial) = (quadratic_observable(V, n-1, n+1)+quadratic_observable(V, n, n+2))/(2*1i);
        end
%         CMI_filename = ['CMI_valuesBeta', num2str(beta_R), 'Mu', num2str(mu_R), 'Dis', num2str(disorder), '.mat'];
%         MI_filename = ['MI_valuesBeta', num2str(beta_R), 'Mu', num2str(mu_R), 'Dis', num2str(disorder), '.mat'];
%         save(CMI_filename, 'CMI_values');
%         save(MI_filename, 'MI_values');
    end
end


toc;
%% Plot
% figure;
% scatter(2*[1:n/2-1], log(abs(MI_values)), 'b');
% xlabel('Size of region B');
% ylabel('Log of Mutual information MI(A, C)');
% 
% figure;
% scatter(2*[1:n/2-1], log(abs(CMI_values)), 'r');
% xlabel('Size of region B');
% ylabel('Log of Conditional Mutual information CMI(A:C|B)');