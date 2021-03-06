function [ current ] = current_analytical(w, gamma, beta_L, mu_L, beta_R, mu_R)
% Implements the analytical formula for the current in a two-site model as
% given in Section D of the transport paper. The current is taken to be I =
% -a_y without w, in order to match the units of the current from Prosen's
% method

r0 = 1/2*(exp(-beta_L*mu_L) + exp(-beta_R*mu_R));
delta_r = 1/2*(exp(-beta_L*mu_L) - exp(-beta_R*mu_R));

current = -2*w/(gamma*(1 + r0)) * delta_r/(4*w^2/gamma^2 + (1+r0)^2 - delta_r^2);
