function [ J ] = compute_J( n )
% Compute the 4nx4n matrix J as defined in Eq. (30)

J = zeros(4*n, 4*n);
for i = 1:2*n
    J(2*i-1, 2*i) = 1;
    J(2*i, 2*i-1) = 1;
end
end

