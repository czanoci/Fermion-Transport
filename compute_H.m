function [ H ] = compute_H( w, n )
% Compute the Hamiltonian for a chain with n sites.
% H is given by the free hopping model
% H = -w * sum(c_x^\dagger c_x+1 + c_x+1^\dagger c_x)

H = zeros(2*n, 2*n);
for x = 1:n-1
    H(2*x-1, 2*x+2) = 1i*w/4;
    H(2*x+2, 2*x-1) = -1i*w/4;
    H(2*x+1, 2*x) = 1i*w/4;
    H(2*x, 2*x+1) = -1i*w/4;
end

end

