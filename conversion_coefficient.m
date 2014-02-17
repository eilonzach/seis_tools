function [CC,EE] = conversion_coefficient(R1,R2,rayp,alpha1,beta1,rho1,alpha2,beta2,rho2)
%[CC,EE] = conversion_coefficient(R1,R2,rayp,alpha1,beta1,rho1,alpha2,beta2,rho2)
%
% CC is the amplitude converseion coefficient
% EE is the energy flux coefficient
%
% rayp is the ray parameter - units of secs/m
% CAREFUL - in spherical earth must divide by radius to get this
% CARFUL with units - velocities and densities must match + think about the
% units of the radius (eg. is asind(rvelocity*rayp0/radius) a realistic angle?
%
% alpha1,beta1,rho1 are P and S velocities and density for the upper medium
% alpha2,beta2,rho2 are P and S velocities and density for the lower medium
%
% R1 and R2 are ray defitions (same as the ones in the path file)
%Change ray definitions used in path file to the ray definitions I used for
% %this script...
% 
% Author C. Eddy  18 May 2013
% Edited Z. Eilon 23 May 2013
% Edited C. Eddy  11 June 2013

if (R1==2)
    r1=2;
elseif (R1==-2)
    r1=1;
elseif (R1==3)
    r1=4;
elseif (R1==-3)
    r1=3;
elseif (R1==1)
    r1=6;
elseif (R1==-1)
    r1=5;
end
if (R2==2)
    r2=2;
elseif (R2==-2)
    r2=1;
elseif (R2==3)
    r2=4;
elseif (R2==-3)
    r2=3;
elseif (R2==1)
    r2=6;
elseif (R2==-1)
    r2=5;
end  

%r1 and r2 are ray definitions for incident and reflected/transmitted wave
%P up = 1
%P down = 2
%SV up = 3
%SV down = 4
%SH up = 5
%SH down = 6   

%%%Calculate angles from ray parameters
i1 = asin(rayp*alpha1);
i2 = asin(rayp*alpha2);
j1 = asin(rayp*beta1);
j2 = asin(rayp*beta2);

%Free surface reflections
%At free surface, assume alpha1=beta1=rho1=0
if (alpha1==0 && beta1==0 && rho1==0)
    eta_alpha = cos(i2)/alpha2;
    eta_beta = cos(j2)/beta2;
if (r1==1 && r2==2) %Upgoing P -> Downgoing P
    CC = (-(1/(beta2^2)-2*rayp^2)^2+4*rayp^2*eta_alpha*eta_beta)/((1/(beta2^2)-2*rayp^2)^2+4*rayp^2*eta_alpha*eta_beta);
    EE = CC;
elseif (r1==1 && r2==4) %Upgoing P -> Downgoing SV
    CC = (4*(alpha2/beta2)*rayp*eta_alpha*(1/(beta2^2)-2*rayp^2))/((1/(beta2^2)-2*rayp^2)^2+4*rayp^2*eta_alpha*eta_beta);
    EE = CC*sqrt((beta2*cos(j2))/(alpha2*cos(i2)));
elseif (r1==3 && r2==4) %Upgoing SV -> Downgoing SV
    CC = ((1/(beta2^2)-2*rayp^2)^2-4*rayp^2*eta_alpha*eta_beta)/((1/(beta2^2)-2*rayp^2)^2+4*rayp^2*eta_alpha*eta_beta);
    EE = CC;
elseif (r1==3 && r2==2) %Upgoing SV -> Downgoing P
    CC = (4*beta2/alpha2*rayp*eta_beta*(1/(beta2^2)-2*rayp^2))/((1/(beta2^2)-2*rayp^2)^2+4*rayp^2*eta_alpha*eta_beta);
    EE = CC*sqrt((alpha2*cos(i2))/(beta2*cos(j2)));
elseif (r1==5 || r2==6) %Upgoing SH -> Downgoing SH
    CC = 1;
    EE = CC;
end

%%%SH waves
%Transmitted 
elseif (r1==5 || r1==6 && rho1~=0) 
if (r1==6 && r2==6) %Downgoing SH -> Downgoing SH
    CC = (2*rho1*beta1*cos(j1))/(rho1*beta1*cos(j1)+rho2*beta2*cos(j2));
    EE = CC*sqrt((rho2*beta2*cos(j2))/(rho1*beta1*cos(j1)));
elseif (r1==5 && r2==5) %Upgoing SH -> Upgoing SH
    CC = (2*rho2*beta2*cos(j2))/(rho1*beta1*cos(j1)+rho2*beta2*cos(j2));
    EE = CC*sqrt((rho1*beta1*cos(j1))/(rho2*beta2*cos(j2)));
%Reflected 
elseif (r1==6 && r2==5) %Downgoing SH -> Upgoing SH
    CC = (rho1*beta1*cos(j1)-rho2*beta2*cos(j2))/(rho1*beta1*cos(j1)+rho2*beta2*cos(j2));
    EE = CC;
elseif (r1==5 && r2==6) %Upgoing SH -> Downgoing SH
    CC = (rho2*beta2*cos(j2)-rho1*beta1*cos(j1))/(rho1*beta1*cos(j1)+rho2*beta2*cos(j2));
    EE = CC;
end

%%P-SV waves

% %Solid-liquid boundary (s above l)
% elseif (beta2==0 && alpha2>0)
%     D = [ cos(i1)   sin(j1)     cos(i2)
%           sin(2*i1)     -(alpha1/beta1)*cos(2*j1)   0
%          -cos(2*j1)     -(beta1/alpha1)*sin(2*j1)   (rho2/rho1)*(alpha2/alpha1)]; 
%     DD = inv(D);
%     v1 = [cos(i1); sin(2*i1); cos(2*j1)];
%     E = [-cos(i1)   -sin(j1)    -cos(i2)
%           (beta1/alpha1)*sin(2*i1)   -cos(2*j1)  0
%           (alpha1/beta1)*cos(2*j1)    sin(2*j1)  -(rho2/rho1)*(alpha2/beta1)];
%     EE = inv(E);
%     v2 = [sin(j1); cos(2*j1); sin(2*j1)];
%     F = [ cos(i2)   cos(i1)     -sin(j1)
%           0         sin(2*i1)   (alpha1/beta1)*cos(2*j1)
%          -1         (rho1/rho2)*(alpha1/alpha2)*cos(2*j1) -(rho1/rho2)*(beta1/alpha2)*sin(2*j1)];
%     FF = inv(F);
%     v3 = [cos(i2); 0; 1];
% if     (r1==2 && r2==1) %Downgoing P -> Upgoing P
%     CC = DD(1,:)*v1;
% elseif (r1==2 && r2==3) %Downgoing P -> Upgoing SV
%     CC = DD(2,:)*v1;
% elseif (r1==2 && r2==2) %Downgoing P -> Downgoing P
%     CC = DD(3,:)*v1;
% elseif (r1==4 && r2==1) %Downgoing SV -> Upgoing P
%     CC = EE(1,:)*v2;
% elseif (r1==4 && r2==3) %Downgoing SV -> Upgoing SV
%     CC = EE(2,:)*v2;
% elseif (r1==4 && r2==2) %Downgoing SV -> Downgoing P
%     CC = EE(3,:)*v2;
% elseif (r1==1 && r2==1) %Upgoing P -> Upgoing P
%     CC = FF(2,:)*v3;
% elseif (r1==1 && r2==3) %Upgoing P -> Upgoing SV
%     CC = FF(3,:)*v3;
% elseif (r1==1 && r2==2) %Upgoing P -> Downgoing P
%     CC = FF(1,:)*v3;
% end
% 
% %Liquid-solid boundary (l above s)
% elseif (beta1==0  && alpha1>0)
%     D = [ cos(i2)   sin(j2)     cos(i1)
%           sin(2*i2)     -(alpha2/beta2)*cos(2*j2)   0
%          -cos(2*j2)     -(beta2/alpha2)*sin(2*j2)   (rho1/rho2)*(alpha1/alpha2)];
%     DD = inv(D);
%     v1 = [cos(i2); sin(2*i2); cos(2*j2)];
%     E = [-cos(i2)   -sin(j2)    -cos(i1)
%           (beta2/alpha2)*sin(2*i2)   -cos(2*j2)  0
%           (alpha2/beta2)*cos(2*j2)    sin(2*j2)  -(rho1/rho2)*(alpha1/beta2)];
%     EE = inv(E);
%     v2 = [sin(j2); cos(2*j2); sin(2*j2)];
%     F = [ cos(i1)   cos(i2)     -sin(j2)
%           0         sin(2*i2)   (alpha2/beta2)*cos(2*j2)
%          -1         (rho2/rho1)*(alpha2/alpha1)*cos(2*j2) -(rho2/rho1)*(beta2/alpha1)*sin(2*j2)];
%     FF = inv(F);
%     v3 = [cos(i1); 0; 1];
% if     (r1==2 && r2==1) %Downgoing P -> Upgoing P
%     CC = FF(1,:)*v3;
% elseif (r1==2 && r2==2) %Downgoing P -> Downgoing P
%     CC = FF(2,:)*v3;
% elseif (r1==2 && r2==4) %Downgoing P -> Downgoing SV
%     CC = FF(3,:)*v3;
% elseif (r1==1 && r2==1) %Upgoing P -> Upgoing P
%     CC = DD(3,:)*v1;
% elseif (r1==1 && r2==2) %Upgoing P -> Downgoing P
%     CC = DD(1,:)*v1;
% elseif (r1==1 && r2==4) %Upgoing P -> Downgoing SV
%     CC = DD(2,:)*v1;
% elseif (r1==3 && r2==1) %Upgoing SV -> Upgoing P
%     CC = EE(3,:)*v2;
% elseif (r1==3 && r2==2) %Upgoing SV -> Downgoing P
%     CC = EE(1,:)*v2;
% elseif (r1==3 && r2==4) %Upgoing SV -> Downgoing SV
%     CC = EE(2,:)*v2;
% end

%Solid-solid boundary
elseif (rho1~=0)
    a = rho2*(1 - 2*beta2.^2*rayp.^2) - rho1*(1 - 2*beta1.^2*rayp.^2);
    b = rho2*(1 - 2*beta2.^2*rayp.^2) + 2*rho1*beta1.^2*rayp.^2;
    c = rho1*(1 - 2*beta1.^2*rayp.^2) + 2*rho2*beta2.^2*rayp.^2;
    d = 2*(rho2*beta2.^2 - rho1*beta1.^2);
    
    E = b*cos(i1)/alpha1 + c*cos(i2)/alpha2;
    F = b*cos(j1)/beta1  + c*cos(j2)/beta2;
    G = a - d*(cos(i1)/alpha1)*(cos(j2)/beta2);
    H = a - d*(cos(i2)/alpha2)*(cos(j1)/beta1);
    D = E*F + G*H*rayp.^2;
    
if (r1==2 && r2==1) %Downgoing P -> Upgoing P
    CC = ((b*cos(i1)/alpha1 - c*cos(i2)/alpha2)*F - (a + d*(cos(i1)/alpha1)*(cos(j2)/beta2))*H*rayp^2)/D;
    EE = CC;
elseif (r1==2 && r2==3) %Downgoing P -> Upgoing SV
    CC = -2*(cos(i1)/alpha1)*(a*b + c*d*(cos(i2)/alpha2)*(cos(j2)/beta2))*rayp*alpha1/(beta1*D);
    EE = CC*sqrt((beta1*cos(j1))/(alpha1*cos(i1)));
elseif (r1==2 && r2==2) %Downgoing P -> Downgoing P
    CC = 2*rho1*cos(i1)/alpha1*F*alpha1/(alpha2*D);
    EE = CC*sqrt((rho2*alpha2*cos(i2))/(rho1*alpha1*cos(i1)));
elseif (r1==2 && r2==4) %Downgoing P -> Downgoing SV
    CC = 2*rho1*cos(i1)/alpha1*H*rayp*alpha1/(beta2*D);
    EE = CC*sqrt((rho2*beta2*cos(j2))/(rho1*alpha1*cos(i1)));
elseif (r1==4 && r2==1) %Downgoing SV -> Upgoing P
    CC = -2*cos(j1)/beta1*(a*b+c*d*cos(i2)/alpha2*cos(j2)/beta2)*rayp*beta1/(alpha1*D);
    EE = CC*sqrt((alpha1*cos(i1))/(beta1*cos(j1)));
elseif (r1==4 && r2==3) %Downgoing SV -> Upgoing SV
    CC = -((b*cos(j1)/beta1-c*cos(j2)/beta2)*E-(a+d*cos(i2)/alpha2*cos(j1)/beta1)*G*rayp^2)/D;
    EE = CC;
elseif (r1==4 && r2==2) %Downgoing SV -> Downgoing P
    CC = -2*rho1*cos(j1)/beta1*G*rayp*beta1/(alpha2*D);
    EE = CC*sqrt((rho2*alpha2*cos(i2))/(rho1*beta1*cos(j1)));
elseif (r1==4 && r2==4) %Downgoing SV -> Downgoing SV
    CC = 2*rho1*cos(j1)/beta1*E*beta1/(beta2*D);
    EE = CC*sqrt((rho2*beta2*cos(j2))/(rho1*beta1*cos(j1)));
elseif (r1==1 && r2==1) %Upgoing P -> Upgoing P
    CC = 2*rho2*cos(i2)/alpha2*F*alpha2/(alpha1*D);
    EE = CC*sqrt((rho1*alpha1*cos(i1))/(rho2*alpha2*cos(i2)));
elseif (r1==1 && r2==3) %Upgoing P -> Upgoing SV
    CC = -2*rho2*cos(i2)/alpha2*G*rayp*alpha2/(beta1*D);
    EE = CC*sqrt((rho1*beta1*cos(j1))/(rho2*alpha2*cos(i2)));
elseif (r1==1 && r2==2) %Upgoing P -> Downgoing P
    CC = -((b*cos(i1)/alpha1-c*cos(i2)/alpha2)*F+(a+d*cos(i2)/alpha2*cos(j1)/beta1)*G*rayp^2)/D;
    EE = CC;
elseif (r1==1 && r2==4) %Upgoing P -> Downgoing SV
    CC = 2*cos(i2)/alpha2*(a*c+b*d*cos(i1)/alpha1*cos(j1)/beta1)*rayp*alpha2/(beta2*D);
    EE = CC*sqrt((beta2*cos(j2))/(alpha2*cos(i2)));
elseif (r1==3 && r2==1) %Upgoing SV -> Upgoing P
    CC = 2*rho2*cos(j2)/beta2*H*rayp*beta2/(alpha1*D);
    EE = CC*sqrt((rho1*alpha1*cos(i1))/(rho2*beta2*cos(j2)));
elseif (r1==3 && r2==3) %Upgoing SV -> Upgoing SV
    CC = 2*rho2*cos(j2)/beta2*E*beta2/(beta1*D);
    EE = CC*sqrt((rho1*beta1*cos(j1))/(rho2*beta2*cos(j2)));
elseif (r1==3 && r2==2) %Upgoing SV -> Downgoing P
    CC = 2*cos(j2)/beta2*(a*c+b*d*cos(i1)/alpha1*cos(j1)/beta1)*rayp*beta2/(alpha2*D);
    EE = CC*sqrt((alpha2*cos(i2))/(beta2*cos(j2)));
elseif (r1==3 && r2==4) %Upgoing SV -> Downgoing SV
    CC = ((b*cos(j1)/beta1-c*cos(j2)/beta2)*E+(a+d*cos(i1)/alpha1*cos(j2)/beta2)*H*rayp^2)/D;
    EE = CC;
end
end