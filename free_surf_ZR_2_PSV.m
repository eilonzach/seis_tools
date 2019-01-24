function [ P,SV ] = free_surf_ZRT_2_PSVSH(Z, R,rayp,alpha,beta,Zconv )
%[ P,SV ] = free_surf_corr( Z, R,rayp,alpha,beta )
%   Function to convert Z and R traces to P,SV
%  R defined as positive away from source
%  Z default is defined as positive upwards; can be switched if Zconv==2
% From Bostock and Rondenay 1999, GJI, we have:
% 
%         (Z->P   R->P    0)
%  invF = (Z->SV  R->SV   0)
%         (  0      0     2)
% 
%         (P->Z   SV->Z   0)
%     F = (P->R   SV->R   0)
%         (  0      0     2)
% 
% where if  u = [Z,R,T]'  and  w = [P,SV,SH]'
% 
%    w = invF*u   --- if w, u are [Nx3] column matrices:   W = U*(invF)'
% or    
%    u = F*w      --- if w, u are [Nx3] column matrices:   U = W*F'

if nargin<6 || isempty(Zconv)
    Zconv=1;
end

if Zconv==1 % Z is defined as positive up
    zsng = 1;
elseif Zconv==2 % Z is defined as positive down
    zsng = -1;
end


rayp = rayp_sdeg2skm(12);
alpha = 6.3;
beta = 3.8;

qa = sqrt(alpha.^-2 - rayp.^2);
qb = sqrt(beta.^-2 - rayp.^2);
DD = 1 - 4*rayp^2*beta^2 + 4*rayp^4*beta^4 + 4*beta^4*rayp^2*qa*qb;

F = zeros(3,3); 
F(3,3)=2;
Fz_p = zsng*2*alpha*qa*(1-2*rayp^2*beta^2)/DD;
Fr_p = 4*alpha*beta^2*rayp*qa*qb/DD;
Fz_sv = zsng*-4*beta^3*rayp*qa*qb/DD;
Fr_sv = 2*beta*qb*(1-2*rayp^2*beta^2)/DD;
F(1:2,1:2) = [Fz_p,Fz_sv;Fr_p,Fr_sv];

invF = inv(F);


end

