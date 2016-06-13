function [ M ] = compute_M( L )
% L is a matrix whose rows are the bath operators represented in the w
% basis. size(L) = (num_baths x 2n)
% M is a (2n x 2n) complex Hermitian matrix given by Eq.(23)

[num_baths, N] = size(L);
% N = 2n

M = zeros(N, N);

for j=1:N
    for k=1:N
        for mu=1:num_baths
           M(j, k) = M(j, k) + L(mu, j)*conj(L(mu, k));
        end
    end
end

end

