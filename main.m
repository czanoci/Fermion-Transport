% Parameters
w = 0.1; % why values w>gamma don't match??
gamma = 3;

t_values = zeros(1, 1000);
c1_values = zeros(1, 1000);
c2_values = zeros(1, 1000);
for i=1:1000
   t = 0.5 + (i-1)*0.1;
   t_values(i) = t;
   c1_values(i) = current_method1(w, gamma, t);
   c2_values(i) = real(current_method2(w, gamma, t));
end
figure;
plot(t_values, c1_values, 'b');
xlabel('Parameter t');
ylabel('Current');
hold on;
plot(t_values, c2_values, 'r');
legend('analytical formula', 'Prosen method');