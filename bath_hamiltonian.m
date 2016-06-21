function [ H ] = bath_hamiltonian( w, bath_size )
% Returns the Hamiltonian for the part of the wire of length bath_size
% written in c_j basis

H = zeros(bath_size, bath_size);

for x=1:bath_size-1
    H(x, x+1) = -w;
    H(x+1, x) = -w;
end


end

