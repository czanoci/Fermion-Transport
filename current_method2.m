function [ current ] = current_method2(w, gamma, beta_L, mu_L, beta_R, mu_R)
% Computes the current for the two-site model using the method described in
% Prosen's paper

%% Define jump operators for the two-site model
uL = sqrt(gamma*exp(-beta_L*mu_L));
uR = sqrt(gamma*exp(-beta_R*mu_R));
vL = sqrt(gamma);
vR = sqrt(gamma);

%% Write the Hamiltonian of the system in the (w1, w2, w3, w4) basis
% See Prosen Eq. (2) 
H = [[0, 0, 0, 1i*w/4];
    [0, 0, -1i*w/4, 0];
    [0, 1i*w/4, 0, 0];
    [-1i*w/4, 0, 0, 0]];

%% Matrix of Lindblad operators
% where each row is a bath operator L_mu written in the (w1, w2, w3, w4) basis 
% L_j's are multiplied 1/sqrt(2) to account for diff. in Lindblad eq. in
% Prosen and transport papers. The eq. in Prosen has an extra factor of 2!
% See Prosen Eq. (3)
L = 1/sqrt(2)*[[0, 0, uR/2, -1i*uR/2];
    [uL/2, -1i*uL/2, 0, 0];
    [0, 0, vR/2, 1i*vR/2];
    [vL/2, 1i*vL/2, 0, 0]];

%% Compute the matrix M defined in Prosen Eq. (23)
M = compute_M(L);

%% Compute the matrix A defined in Prosen Eq. (27)
A = compute_A(H, M);

% Make sure matrix A is anti-symmetric
if A ~= -A.'
    disp('Matrix A not anti-symmetric');
end

%% Compute the current using Theorem 3 from Prosen's paper
current = compute_current(A);

