function [ current ] = current_method2(w, gamma, beta_L, mu_L, beta_R, mu_R)

uL = sqrt(gamma*exp(-beta_L*mu_L));
uR = sqrt(gamma*exp(-beta_R*mu_R));
vL = sqrt(gamma);
vR = sqrt(gamma);

% Hamiltonian
H = [[0, 0, 0, 1i*w/4];
    [0, 0, -1i*w/4, 0];
    [0, 1i*w/4, 0, 0];
    [-1i*w/4, 0, 0, 0]];

% Matrix of Lindblad operators ADDED 1/sqrt(2) to account for diff in
% lindblad eq.
L = 1/sqrt(2)*[[0, 0, uR/2, -1i*uR/2];
    [uL/2, -1i*uL/2, 0, 0];
    [0, 0, vR/2, 1i*vR/2];
    [vL/2, 1i*vL/2, 0, 0]];

M = compute_M(L);
A = compute_A(H, M);

if A ~= -A.'
    disp('Matrix A not anti-symmetric');
end

[current, V, eig] = compute_current(A);

