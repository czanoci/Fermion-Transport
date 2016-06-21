function [ L ] = compute_L( w, gamma, beta_L, mu_L, beta_R, mu_R, nL, nW, nR )
% Computes the matrix L of Lindblad operators
% each row represents a bath operator written in w_j basis

hL = bath_hamiltonian( w, nL );
hR = bath_hamiltonian( w, nR );

[eigenvectors_L, eigenvalues_L] = eig(hL);
eigenvalues_L = diag(eigenvalues_L);
[eigenvectors_R, eigenvalues_R] = eig(hR);
eigenvalues_R = diag(eigenvalues_R);

L = zeros(2*(nL+nR), 2*(nL+nW+nR));

for row = 1:nR
    uR = sqrt(gamma*exp(beta_R*(eigenvalues_R(row) - mu_R)));
    vR = sqrt(gamma);
    psi = eigenvectors_R(:, row);
    for col = 1:nR
        L(row, 2*nL + 2*nW + 2*col - 1) = 1/2*uR*psi(col);
        L(row, 2*nL + 2*nW + 2*col) = -1i/2*uR*psi(col);
        L(nR + row, 2*nL + 2*nW + 2*col - 1) = 1/2*vR*psi(col);
        L(nR + row, 2*nL + 2*nW + 2*col) = 1i/2*vR*psi(col);
    end
end

for row = 1:nL
    uL = sqrt(gamma*exp(beta_L*(eigenvalues_L(row) - mu_L)));
    vL = sqrt(gamma);
    psi = eigenvectors_L(:, row);
    for col = 1:nL
        L(2*nR + row, 2*col - 1) = 1/2*uL*psi(col);
        L(2*nR + row, 2*col) = -1i/2*uL*psi(col);
        L(2*nR + nL + row, 2*col - 1) = 1/2*vL*psi(col);
        L(2*nR + nL + row, 2*col) = 1i/2*vL*psi(col);
    end
end

L = 1/sqrt(2)*L;
