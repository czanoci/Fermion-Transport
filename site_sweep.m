%% Parameters
w = 1; 
ww = 1;
gamma = 1;
beta_L = 5;
beta_R = 5;
mu_L = -0.01;
mu_R = 0.01;
nL = 25;
nR = 25;
nW = 50;

tic;
%% Evaluation
current_values = zeros(1, nL+nR+nW-1);
site_values = zeros(1, nL+nR+nW-1);

for site=1:nL+nR+nW-1
    current = compute_current_multisite(w, ww, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR, site, 0);
    disp(site);
    disp(current);
    if site == nL || site == nL+nW
        current_values(site) = ww*real(current);
    else 
        current_values(site) = w*real(current);
    end
    site_values(site) = site;
end

toc;
%% Plot
figure;
plot(site_values, current_values, 'b');
xlabel('Measurement site');
ylabel('Current');