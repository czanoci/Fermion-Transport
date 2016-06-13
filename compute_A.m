function [ A ] = compute_A( H, M)
% H - (2n x 2n) Hamiltonian matrix in w basis
% M - (2n x 2n) Hermitian matrix as outputed by compute_M
% A - (4n x 4n) complex antisummetric matrix computed according to Eq. (27)

% N = 2n
[N, ~] = size(H);
A = zeros(2*N, 2*N);

for j=1:N
    for k=1:N
        A(2*j-1, 2*k-1) = -2*1i*H(j, k) - M(j, k) + M(k, j);
        A(2*j-1, 2*k) = 2*1i*M(k, j);
        A(2*j, 2*k-1) = -2*1i*M(j, k);
        A(2*j, 2*k) = -2*1i*H(j, k) + M(j, k) - M(k, j);
    end
end

end

