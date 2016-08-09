function [ H_aa ] = aubry_andre_model( V, nL, nW, nR, Q, theta )
% Add the Aubry-Andre term to the Hamiltonian for a chain divided into:
% 1) Left bath with length nL
% 2) The actual wire with length nW
% 3) Right bath with length nR
% The cosine potential applies only to the nW part of the wire and is given by
% H_dis = sum( V cos(Q*n+theta) c_r^\dagger c_r)
% where Q is irrational
% NOTE: The Hamiltonian is written in w_j, not c_j basis.

n = nL + nW + nR;
H_aa = zeros(2*n, 2*n);
for r = nL+1:nL+nW
    V_r = V*cos(Q*r+theta);
    %H_aa(2*r-1, 2*r-1) = V_r/4;
    H_aa(2*r, 2*r-1) = 1i*V_r/4;
    H_aa(2*r-1, 2*r) = -1i*V_r/4;
    %H_aa(2*r, 2*r) = V_r/4;
end

end

