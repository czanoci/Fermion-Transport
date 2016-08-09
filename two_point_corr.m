function [ corr ] = two_point_corr( r, s, V )
% Computes the expectation value <c_r^\dagger c_s> for lattice sites r and
% s

corr = (quadratic_observable(V, 2*r-1, 2*s-1) + quadratic_observable(V, 2*r, 2*s) + 1i*quadratic_observable(V, 2*r, 2*s-1) - 1i*quadratic_observable(V, 2*r-1, 2*s))/4;

end

