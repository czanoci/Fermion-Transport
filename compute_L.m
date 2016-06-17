function [ L ] = compute_L( uR, uL, vR, vL, n )
% Computes the matrix L of Lindblad operators

L = zeros(4, 2*n);

L(1, 2*n-1) = uR/(2*sqrt(2));
L(1, 2*n) = -1i*uR/(2*sqrt(2));

L(2, 1) = uL/(2*sqrt(2));
L(2, 2) = -1i*uL/(2*sqrt(2));

L(3, 2*n-1) = vR/(2*sqrt(2));
L(3, 2*n) = 1i*vR/(2*sqrt(2));

L(4, 1) = vL/(2*sqrt(2));
L(4, 2) = 1i*vL/(2*sqrt(2));

end

