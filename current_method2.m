function [ current ] = current_method2( w, gamma, t)

uL = sqrt(gamma*t);
uR = sqrt(gamma/t);

% Hamiltonian
H = [[0, 0, 0, 1i*w/4];
    [0, 0, -1i*w/4, 0];
    [0, 1i*w/4, 0, 0];
    [-1i*w/4, 0, 0, 0]];

% Matrix of Lindblad operators
L = [[0, 0, uR/2, -1i*uR/2];
    [uL/2, -1i*uL/2, 0, 0];
    [0, 0, sqrt(gamma)/2, 1i*sqrt(gamma)/2];
    [sqrt(gamma)/2, 1i*sqrt(gamma)/2, 0, 0]];

M = compute_M(L);
A = compute_A(H, M);

if A ~= -A.'
    disp('Matrix A not anti-symmetric');
end

[current, V, eig] = compute_current(A);




% [n, ~] = size(V);
% N = zeros(n, n);
% for i=1:n
%     for j=1:n
%         N(i, j) = V(i, :)*V(j, :).';
%     end
% end
% disp(N);






