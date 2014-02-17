function [dtN,dtE,dtNE,acorN,acorE,acorNE,stdN,stdE,stdNE] = SEN_xcortimes(datNE,dt,pretime,lagmax,plotopt)
% [dtN,dtE,dtNE,acorN,acorE,acorNE,stdN,stdE,stdNE] = SEN_xcortimes(datNE,dt,pretime,lagmax,plotopt)
% Augmented version of Vandecar and Crosson for shear waves on two
% components, assuming the two components are the null directions of the
% anisotropy, and hence the two shear pulses will be identically shaped,
% but phase-shifted by a delay time

%% FOR TESTING
% % nstas = 3;
% % dtNE = random('unif',0.3,0.7,2,nstas);
% % datNE = zeros(100,nstas,2);
% % for is = 1:nstas
% %     for ic = 1:2
% %     datNE(:,is,ic) = synthtrace(50,5,10,0.5,'gauss',dtNE(ic,is));
% %     end
% % end
% load('trace_example.mat')
% tt = [0:0.02:32];
% figure(1)
% subplot(2,1,1)
% plot(tt,datNE(:,:,1));
% subplot(2,1,2)
% plot(tt,datNE(:,:,2));
% 
% pretime = 8;
% lagmax = 30;
% dt = 1/50;
% plotopt = 7;
% nargin=5; fprintf('remember to delete nargin setting\n')
%% END TESTING


if nargin<2
    dt=1; % default is 1 sps
end
if nargin<3
    pretime=length(datNE(:,1,1))*dt/2; %default is pick is halfway 
end
if nargin<4
    lagmax = length(datNE(:,1,1))*dt; % default is no lagmax
end
if nargin<5
    plotopt=0; % default is no plot
end


nstas = size(datNE,2);
npts = size(datNE,1);

N = handshake(nstas)*2 + nstas + 1;
M = nstas*2;
G = zeros(N,M);
d = zeros(N,1);
acor = zeros(N,1);

lmax=round(lagmax./dt);

kk = 1;
for ii = 1:(nstas-1)
    for jj = (ii+1):nstas
        for cc = 1:2
        dz = handshake(nstas)*(cc-1);
        dx = nstas*(cc-1);
        [c,lags]=xcorr(datNE(:,ii,cc),datNE(:,jj,cc),lmax,'biased');
        pk=find(c==max(c));
        d(dz + kk)  = lags(pk(1));
        acor(dz + kk) = max(c)/(std(datNE(:,ii,cc))*std(datNE(:,jj,cc)));
        G(dz + kk,dx + ii) =  1;
        G(dz + kk,dx + jj) = -1;
        end
        kk=kk+1;       
    end
end

for ii = 1:nstas
    dz = handshake(nstas)*2;
    dx = nstas;
    [c,lags]=xcorr(datNE(:,ii,1),datNE(:,ii,2),lmax,'biased');
    pk=find(c==max(c));
    d(dz + ii)  = lags(pk(1));
    acor(dz + ii) = max(c)/(std(datNE(:,ii,1))*std(datNE(:,ii,2)));
    G(dz + ii,ii)      =  1;
    G(dz + ii,dx + ii) = -1;
end

%% optional - G is not singular as is...
% decide that mean of N arrivals is zero
 G(end,1:nstas) = 1;
 acor(end) = 1;


%% Solve for m
Wd = diag(acor);
% Wd = eye(size(Wd));
m = (G'*Wd*G)\G'*Wd*d;
resid = d - G*m;
mN = m(1:nstas);
mE = m(nstas+(1:nstas));

%% Tdiff results
dtN = mN*dt + pretime;
dtE = mE*dt + pretime;
dtNE = dtN-dtE;

[dtN dtE dtNE]

%% find acor
stakN = zeros(npts,1);
stakE = zeros(npts,1);
tstak = [0:npts-1]';
for is=1:nstas
    stakN=stakN+interp1(tstak-mN(is),datNE(:,is,1),tstak);
    stakE=stakE+interp1(tstak-mE(is),datNE(:,is,2),tstak);
end
stakN=stakN./nstas;
stakE=stakE./nstas;

% AcorN
jmx=npts-ceil(max(mN));
jmn=ceil(abs(min(mN)))+1;
ida=find(isfinite(stakN) & (1:npts)'<jmx & (1:npts)'>jmn );
acorN = zeros(nstas,1);
for is=1:nstas
    acorN(is)=sum(stakN(ida).*datNE(ida+round(mN(is)),is,1))./(std(stakN(ida),1).*std(datNE(ida+round(mN(is)),is,1),1))./length(ida);
end

% AcorE
jmx=npts-ceil(max(mE));
jmn=ceil(abs(min(mE)))+1;
ida=find(isfinite(stakE) & (1:npts)'<jmx & (1:npts)'>jmn );
acorE = zeros(nstas,1);
for is=1:nstas
    acorE(is)=sum(stakE(ida).*datNE(ida+round(mE(is)),is,2))./(std(stakE(ida),1).*std(datNE(ida+round(mE(is)),is,2),1))./length(ida);
end

acorNE = acor(handshake(nstas)*2 + (1:nstas));

%% Plot
if (plotopt>0) 
    figure(plotopt), clf
tt0=dt.*(0:(npts-1))-pretime;

subplot(223);  plot((tstak*dt)-2*pretime,stakN,'r','Linewidth',5)
subplot(224);  plot((tstak*dt)-2*pretime,stakE,'r','Linewidth',5)


for is=1:nstas
    subplot(221); hold on
    plot(tt0,datNE(:,is,1),'Linewidth',1.5);
    title('Pre-alignment North')
    
    subplot(223); hold on
    plot(tt0-dtN(is),datNE(:,is,1),'Linewidth',1.5);
    title('Post-alignment North')
    
    subplot(222); hold on
    plot(tt0,datNE(:,is,2),'Linewidth',1.5);
    title('Pre-alignment East')
    
    subplot(224); hold on
    plot(tt0-dtE(is),datNE(:,is,2),'Linewidth',1.5);
    title('Post-alignment East')
end

end % on plotopt

end
