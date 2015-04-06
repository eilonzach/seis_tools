function [phi, tlag,ergrid,jmin, rslt]=splitga_sc(dan,dae,dt)
%  function [phi, tlag,ergrid,jmin, rslt]=splitga_sc(datn,date,dt)
%  Calculate splitting parameters from shear wave, following Silver+Chan
%  1991 JGR
%   INPUT:  
%      dan, dae  = N, E component seismograms, already windowed on SKS
%                    assumed to be time-aligned
%      dt = samplerate
%  OUTPUT:
%   phi = array of N  angles tested  ( -90 t  90), degrees
%   tlag = array of M delay times tested ( 0 to 8), sec
%   ergrid = M x N grid of smallest eigenvalues, this is minimized
%   jmin(2) = index of solution; best fit is ergrid(jmin(1),jmin(2))
%   rslt(2) = [tlag(jmin(1)), phi(jmin(2))  the solution
%               

phi=-90:90;
tlag=0:dt:5.0;
N = length(phi);
M = length(tlag);
ergrid=zeros(M, N);
rpd=pi/180;

% Make cross-correlation matrix
len=min([length(dan) length(dae)]);
wdow=tukeywin(len,0.1);
dn=reshape(dan(1:len),len,1).*wdow;
de=reshape(dae(1:len),len,1).*wdow;
dmat=[dn, de];
[cc,lags]=xcorr(dmat);
i0=find(lags==min(abs(lags)));
lags=lags.*dt;

cc0=reshape(cc(i0,:),2,2); % zero-lag cross-correlation; 

ermin=1.e30;
jmin=[0,0];
rslt=[0,0];
for ii=1:N
  ph=(phi(ii)+90.)*rpd;     % add +90. to fix some confusion about axes 10/09
  cp=cos(ph);
  sp=sin(ph);
  rmat=[cp, sp; -sp,cp];
  cbar=(rmat)*cc0*(rmat');   % this has 0-lag for diagonals
  for jj=1:M
    tt=-tlag(jj);
    k=find(abs(lags-tt)==min(abs(lags-tt)));
    cct=reshape(cc(k,:),2,2);
    cb21=rmat(2,1)*(rmat(1,1)*cct(1,1)+rmat(1,2)*cct(2,1)) + rmat(2,2)* ...
	 (rmat(1,1)*cct(1,2)+rmat(1,2)*cct(2,2));
    cbar(2,1)=cb21;
    cbar(1,2)=cb21;
    b=-cbar(1,1)-cbar(2,2);
    c=cbar(1,1)*cbar(2,2)-cbar(1,2)*cbar(2,1);
    radic=(b*b-4*c);
    if (radic<0.) 
      disp('CPLX Eig for dt %.2f phi %.0f',-tt,phi(ii))
    end
    sqrad=sqrt(radic);
    %    d=eig(cbar);
    %lam1=min(abs(d));
    lam1=0.5*(-b-sqrad);
    ergrid(jj,ii)=lam1;
    l1=abs(lam1);
    l2=abs(0.5*(-b+sqrad));
    if (lam1 < ermin & l1 < l2)
      ermin=l1;
      jmin=[jj,ii];
      rslt=[tlag(jj), phi(ii)];
    elseif (l2<ermin & l2 < l1)
      ermin=l2;
      jmin=[jj, ii];
      rslt=[tlag(jj), phi(ii)];
    end
  end
end
ergrid = ergrid ./ ermin;


    