function [ E ] = compute_energy( V, w, nL, nW )
%Compute the energy of the middle part of the wire, i.e. the expected value
%of the hamiltonian

E = 0;

for x = nL+1:nL+nW-1
    disp(E);
    E = E + quadratic_observable( V, 2*x-1, 2*x+2 ) + quadratic_observable( V, 2*x+1, 2*x );
end

E = 1i*w/2*E;

