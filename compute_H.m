function [ H ] = compute_H( w, ww, nL, nW, nR )
% Compute the Hamiltonian for a chain divided into:
% 1) Left bath with length nL
% 2) The actual wire with length nW
% 3) Right bath with length nR
% Each of the 3 parts has a decoupled Hamiltonian
% H = -w * sum(c_x^\dagger c_x+1 + c_x+1^\dagger c_x)
% The different regions are coupled by -ww * (c_nL^\dagger c_nL+1 + c_nL+1^\dagger c_nL)
% NOTE: The Hamiltonian is written in w_j, not c_j basis.

n = nL + nW + nR;
H = zeros(2*n, 2*n);
for x = 1:n-1
    if x == nL || x == nL + nW
        H(2*x-1, 2*x+2) = 1i*ww/4;
        H(2*x+2, 2*x-1) = -1i*ww/4;
        H(2*x+1, 2*x) = 1i*ww/4;
        H(2*x, 2*x+1) = -1i*ww/4;
    else
        H(2*x-1, 2*x+2) = 1i*w/4;
        H(2*x+2, 2*x-1) = -1i*w/4;
        H(2*x+1, 2*x) = 1i*w/4;
        H(2*x, 2*x+1) = -1i*w/4;
    end
end

end

