function [ wjwk ] = quadratic_observable( V, j, k )
% Computes the expected value of quadratic observable <wjwk> using Theorem
% 3 from Prosen's paper (see Eq. (47)).

N = size(V, 1);
n = N/4;

% if j == k
%     wjwk = 1;
% else
%     wjwk = 0;
% end
% 
% for m = 1:2*n
%    wjwk = wjwk + 1/2 * (V(2*m, 2*j-1)*V(2*m-1, 2*k-1) - V(2*m, 2*j)*V(2*m-1, 2*k) - 1i*V(2*m, 2*j)*V(2*m-1, 2*k-1) - 1i*V(2*m, 2*j-1)*V(2*m-1, 2*k));
% end

wjwk = 0;
for m = 1:2*n
    wjwk = wjwk + 2*V(2*m, 2*j-1)*V(2*m-1, 2*k-1);
end
end

