function [ nu ] = poisson_from_Vp_Vs( Vp,Vs )
%[ nu ] = poisson_from_Vp_Vs( Vp,Vs )
%   simple function to compute the effective Poisson's ratio from the
%   values of Vp and Vs

% there should be a factor of rho in the following equations but it cancels
% so we ignore it.
M = Vp.^2;
G = Vs.^2;

nu = (M - 2*G)./(2*M - 2*G); 


end

