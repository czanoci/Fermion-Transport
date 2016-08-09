function [ H_dis ] = add_disorder( V0, nL, nW, nR )
% Add the disorder term to the Hamiltonian for a chain divided into:
% 1) Left bath with length nL
% 2) The actual wire with length nW
% 3) Right bath with length nR
% The disorder applies only to the nW part of the wire and is given by
% H_dis = sum( V_r c_r^\dagger c_r)
% where V_r is uniformly distributed between (-V0, +V0)
% NOTE: The Hamiltonian is written in w_j, not c_j basis.
%rng(42);
n = nL + nW + nR;
H_dis = zeros(2*n, 2*n);
for r = nL+1:nL+nW
    V_r = -V0+2*V0*rand;
    %H_dis(2*r-1, 2*r-1) = V_r/4;
    H_dis(2*r, 2*r-1) = 1i*V_r/4;
    H_dis(2*r-1, 2*r) = -1i*V_r/4;
    %H_dis(2*r, 2*r) = V_r/4;
end

end

