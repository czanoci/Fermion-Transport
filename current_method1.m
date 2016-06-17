function [ current ] = current_method1(w, gamma, beta_L, mu_L, beta_R, mu_R)

r0 = 1/2*(exp(-beta_L*mu_L) + exp(-beta_R*mu_R));
delta_r = 1/2*(exp(-beta_L*mu_L) - exp(-beta_R*mu_R));

current = -2*w/(gamma*(1 + r0)) * delta_r/(4*w^2/gamma^2 + (1+r0)^2 - delta_r^2);
