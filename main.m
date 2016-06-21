%% Parameters
w = 10; 
gamma = 1;
beta_L = 2;
beta_R = 1;

%% Evaluation
t_values = zeros(1, 1000);
c1_values = zeros(1, 1000);
c2_values = zeros(1, 1000);
for i=1:1000
   t = 0.5 + (i-1)*0.1;
   t_values(i) = t;
   c1_values(i) = current_analytical(w, gamma, beta_L, -log(t)-log(gamma), beta_R, log(t)-log(gamma));
   c2_values(i) = real(current_prosen(w, gamma, beta_L, -log(t)-log(gamma), beta_R, log(t)-log(gamma)));
end

%% Plot
figure;
plot(t_values, c1_values, 'b');
xlabel('Parameter t');
ylabel('Current');
hold on;
plot(t_values, c2_values, 'r');
legend('analytical formula', 'Prosen method');