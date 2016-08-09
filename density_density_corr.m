function [ corr ] = density_density_corr( r, s, V )
% Computes the density-density correlation function <c_r^\dagger c_r c_s^\dagger c_s> 
% for lattice sites r and s

if r == s
    corr = two_point_corr( r, s, V );
else
    corr = 0;
end

corr = corr + two_point_corr( r, r, V )*two_point_corr( s, s, V ) - two_point_corr( r, s, V )*two_point_corr( s, r, V );

end

