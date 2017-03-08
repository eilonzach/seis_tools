function p=PropMtxRay(omega,kx,h,rho,alpha,beta) 

% Calcualting Propagator Matrix for Rayleigh wave (P-SV wave) for any sinlge
% layer
% h: thickness
% alpha:vp
% beta: vs 
% P: 4*4 matrix
%
% p=zeros(4,4);

if (kx.^2 - (omega/alpha)^2 > 0) 
    gamma = sqrt( kx.^2 - (omega/alpha)^2);
else 
    gamma = -sqrt((omega/alpha)^2-kx.^2 ) * 1i;    
end

if (sqrt(kx.^2 - (omega/beta)^2) > 0)
    eta = sqrt(kx.^2 - (omega/beta)^2);
else
    eta = -sqrt((omega/beta)^2-kx.^2 ) * 1i;    
    
end

mu=rho*beta^2;  % shear modulus 

% only do these once for speed
w2 = omega^2;
e2 = eta^2;
w2r = w2*rho;
kx2 = kx^2;
kx2_p_e2 = kx2+e2;
shgh = sinh(gamma*h);
sheh = sinh(eta*h);
shgh_2 = sinh(gamma*h/2); shgh_2shgh_2 = shgh_2^2;
sheh_2 = sinh(eta*h/2); sheh_2sheh_2 = sheh_2^2;
m_b_w2r = mu/w2r;
m2_b_w2r = mu^2/w2r;

% Construct propagator matrix for a single layer
% reference, Aki and Richards Quantiative Seismology 2nd Ed. pp 273
p(1,1)=1 + 2*m_b_w2r*(2*kx2*shgh_2shgh_2 - kx2_p_e2*sheh_2sheh_2);
p(1,2)=kx*m_b_w2r*(kx2_p_e2/gamma*shgh - 2*eta*sheh);
p(1,3)=1/w2r*(kx2/gamma*shgh - eta*sheh);
p(1,4)=2*kx/w2r*(shgh_2shgh_2 - sheh_2sheh_2);
p(2,1)=kx*m_b_w2r*(kx2_p_e2/eta*sheh - 2*gamma*shgh);
p(2,2)=1 + 2*m_b_w2r*(2*kx2*sheh_2sheh_2 - kx2_p_e2*shgh_2shgh_2);
p(2,3)=-p(1,4);
p(2,4)=1/w2r*(kx2/eta*sheh - gamma*shgh);
p(3,1)=m2_b_w2r*(4*kx2*gamma*shgh - kx2_p_e2^2/eta*sheh);
p(3,2)=2*mu^2*kx2_p_e2*p(1,4);
p(3,3)=p(1,1);
p(3,4)=-p(2,1);
p(4,1)=-p(3,2);
p(4,2)=m2_b_w2r*(4*kx2*eta*sheh - kx2_p_e2^2/gamma*shgh);
p(4,3)=-p(1,2);
p(4,4)=p(2,2);
return
end
