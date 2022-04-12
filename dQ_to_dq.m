function [ dq ] = dQ_to_dq( dQ )
% [ dq ] = dQ_to_dq( dQ )
% or
% [ dQ ] = dQ_to_dq( dq )
% 
% where dQ and dq are fractional changes in Q or q respectively. 


dq = -dQ./(1+dQ);

end