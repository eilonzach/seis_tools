function resp=freesurf_resp(rp, alph, bet)
% resp=freesurf_resp(rp, alph, bet)
% freesurface correction factor for incident plane wave
% output:  resp(iph,comp)  where 
%       iph=1,2,3 for P,SV, SH of incident wave
%       icomp = 1,2, 3 for Z, R, T,  Z +ve up, R +ve in dir'n of wave propagation
%	subroutine freesurf(rp,alph,bet,iph,icomp,resp)
%	input:	rp = ray parameter (sin(i)/v)
%		alph = receiver p-velocity
%		bet = receiver s-velocity
%  Only works for SCALARS
%
%   G. Abers 2002 

resp=zeros(3,3);
resp(3,3) = 2.0;

ai = asin(alph.*rp);
aj = asin(bet.*rp);
ci = cos(ai);
cj = cos(aj);
c2j = cos(2.*aj);
s2i = sin(2.*ai);
s2j = sin(2.*aj);
rat = bet./alph;
dd = c2j.*c2j + rat.*rat.*s2i.*s2j;     % 7/02:  THIS IS FIXED; was wrong in all older versions
resp(1,1) = 2.*ci.*c2j./dd;
resp(1,2) = 2.*ci.*s2j./dd;
resp(2,2) = 2.*cj.*c2j./dd;
resp(2,1) = -2.*rat.*ci.*s2j./dd;     % 7/02:  FIXED sign convention to be consistent w/ iph=1

return
