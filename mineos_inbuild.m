% Script to make a 1D earth profile that can go into the normal mode
% summation program MINEOS.
% by default the model file is in tabular form
% the model file will appear in this dir as "NAME.model"

name = 'pprem'; % name of model
odir = '~/Documents/MATLAB/Classes/Adv.'' Seismp''/';

ifanis = 1; % 1 if there is some (radial) anisotropy, 0 otherwise
tref = 1.;  % reference period (recommend 1) for physical dispersion

N = 300; % no. of nodes
Ric  = 6371e3 - 5150e3; % Radius of inner/outer core boundary
Roc  = 6371e3 - 2891e3; % Radius of outer core/mantle boundary
Rlm1 = 3630e3;
Rlm2 = 5600e3;
Rtzb = 6371e3 - 670e3; % Radius of base of transition zone
Rtzt = 6371e3 - 440e3; % Radius of top of transition zone
Rlvz = 6371e3 - 220e3; % Radius of base of low-velocity zone
Rlab = 6371e3 - 125e3; % Radius of lith/asth bounary
Rmoh = 6371e3 - 35e3; % Radius of moho
Ruc  = 6371e3 - 20e3; % Radius of upper/lower crust boundary
Rsfl = 6371e3; % Radius of seafloor or land/air
Rsfc = 6371e3; % Radius of sea/air or land/air

Rs = [0 Ric Roc Rlm1 Rlm2 Rtzb Rtzt Rlvz Rlab Rmoh Ruc Rsfl Rsfc]';
Nlay = length(Rs) - 1;

% Coefficients of models - each line is a different layer, and columns are
% coefficients of indices 0, 1, 2, 3 respectively

Dc =   [13.0885      0       -8.8381      0
        12.5815     -1.2638  -3.6426     -5.5281
         7.9565     -6.4761   5.5283     -3.0807
         7.9565     -6.4761   5.5283     -3.0807
         7.9565     -6.4761   5.5283     -3.0807
        10.6449     -7.4347   0           0
         7.1089     -3.8045   0           0
         2.6910     +0.6924   0           0
         2.6910     +0.6924   0           0
         2.9000      0        0           0
         2.6000      0        0           0
         1.0200      0        0           0      ];

Vpvc = [11.2622       0       -6.3640      0
        11.0487      -4.0362   4.8023    -13.5732
        15.3891      -5.3181   5.5242     -2.5514
        24.9520     -40.4673  51.4832    -26.6419
        29.2766     -23.6027   5.5242     -2.5514
        34.1736     -26.7171   0           0
        20.3926     -12.2569   0           0
         0.8317       7.2180   0           0
         0.8317       7.2180   0           0
         6.8000       0        0           0
         5.8000       0        0           0
         1.4500       0        0           0      ];
       
 
Vphc = [11.2622       0       -6.3640      0
        11.0487      -4.0362   4.8023    -13.5732
        15.3891      -5.3181   5.5242     -2.5514
        24.9520     -40.4673  51.4832    -26.6419
        29.2766     -23.6027   5.5242     -2.5514
        34.1736     -26.7171   0           0
        20.3926     -12.2569   0           0
         3.5908       4.6172   0           0
         3.5908       4.6172   0           0
         6.8000       0        0           0
         5.8000       0        0           0
         1.4500       0        0           0      ];   
    
Vsvc = [ 3.6678       0       -4.4475      0
         0        	  0        0           0
         6.9254       1.4672  -2.0834      0.9783
        11.1671     -13.7818  17.4575     -9.2777
        22.3459     -17.2473  -2.0834      0.9783
        19.0335     -15.0456   0           0
         8.9496      -4.4597   0           0
         5.8582      -1.4678   0           0
         5.8582      -1.4678   0           0
         3.9000       0        0           0
         3.2000       0        0           0
         0            0        0           0      ];
     
Vshc = [ 3.6678       0       -4.4475      0
         0        	  0        0           0
         6.9254       1.4672  -2.0834      0.9783
        11.1671     -13.7818  17.4575     -9.2777
        22.3459     -17.2473  -2.0834      0.9783
        19.0335     -15.0456   0           0
         8.9496      -4.4597   0           0
        -1.0839       5.7176   0           0
        -1.0839       5.7176   0           0
         3.9000       0        0           0
         3.2000       0        0           0
         0            0        0           0      ]; 

Qmc = [ 84.6;
        9999999
        312
        312
        312
        143
        143
        80
        600
        600
        600
        9999999 ]; 

Qkc = [ 1327.7;
        57823
        57823
        57823
        57823
        57823
        57823
        57823
        57823
        57823
        57823
        57823 ];   

Etac = [1        0
        1        0
        1        0
        1        0
        1        0
        1        0
        1        0
        3.3687  -2.4778
        3.3687  -2.4778
        1        0
        1        0
        1        0 ];

%     PROBLEM - have to make sure each layer has at least two knots - i.e. round up if below two 
    

% Each layer has a bottom and a top knot.  No. of nodes in layer is therefore
% Ni = 1 + round((Rtop-Rbottom)/Rtotal)*(N-Nlay)
% and spacing within layer is (use interp) (Rtop-Rbottom)/(Ni-1)
Nk = zeros(Nlay,1);
R = zeros(Nlay,1); % radius
D = zeros(Nlay,1); % density
Vpv = zeros(Nlay,1); 
Vph = zeros(Nlay,1);
Vsv = zeros(Nlay,1);
Vsh = zeros(Nlay,1);
Qm = zeros(Nlay,1);
Qk = zeros(Nlay,1);
Eta = zeros(Nlay,1);

ibo = 1;
for kk = 1:Nlay
Nk(kk) = 2 + round(((Rs(kk+1)-Rs(kk))/Rs(end))*(N - 2*Nlay));
ito = sum(Nk(1:kk));

R(ibo:ito) = linspace(Rs(kk),Rs(kk+1),Nk(kk))';
X = R/Rs(end);

D(ibo:ito)   = Dc(kk,1)  *ones(Nk(kk),1) + Dc(kk,2)  *X(ibo:ito) + Dc(kk,3)  *X(ibo:ito).^2 + Dc(kk,4)  *X(ibo:ito).^3;
Vpv(ibo:ito) = Vpvc(kk,1)*ones(Nk(kk),1) + Vpvc(kk,2)*X(ibo:ito) + Vpvc(kk,3)*X(ibo:ito).^2 + Vpvc(kk,4)*X(ibo:ito).^3;
Vsv(ibo:ito) = Vsvc(kk,1)*ones(Nk(kk),1) + Vsvc(kk,2)*X(ibo:ito) + Vsvc(kk,3)*X(ibo:ito).^2 + Vsvc(kk,4)*X(ibo:ito).^3;
Vph(ibo:ito) = Vphc(kk,1)*ones(Nk(kk),1) + Vphc(kk,2)*X(ibo:ito) + Vphc(kk,3)*X(ibo:ito).^2 + Vphc(kk,4)*X(ibo:ito).^3;
Vsh(ibo:ito) = Vshc(kk,1)*ones(Nk(kk),1) + Vshc(kk,2)*X(ibo:ito) + Vshc(kk,3)*X(ibo:ito).^2 + Vshc(kk,4)*X(ibo:ito).^3;
Qm(ibo:ito) = Qmc(kk,1)*ones(Nk(kk),1);
Qk(ibo:ito) = Qkc(kk,1)*ones(Nk(kk),1);
Eta(ibo:ito) = Etac(kk,1)*ones(Nk(kk),1) + Etac(kk,2)*X(ibo:ito);

ibo = 1 + ito;
end

%% Plot up model
figure(1); clf
hold on
plot((Rs(end)-R)/1000,D,'k')
plot((Rs(end)-R)/1000,Vpv,'b')
plot((Rs(end)-R)/1000,Vph,'--b')
plot((Rs(end)-R)/1000,Vsv,'r')
plot((Rs(end)-R)/1000,Vsh,'--r')
xlim([0 Rs(end)/1000])
xlabel('Radius (km)')
ylabel('Density (kg/m3)')
hold off

figure(2); clf
hold on
plot((Rs(end)-R)/1000,D,'k')
plot((Rs(end)-R)/1000,Vpv,'b')
plot((Rs(end)-R)/1000,Vph,'--b')
plot((Rs(end)-R)/1000,Vsv,'r')
plot((Rs(end)-R)/1000,Vsh,'--r')
plot((Rs(end)-R)/1000,Eta,'g')
xlim([0 1e6])
xlabel('Radius (km)')
ylabel('Density (kg/m3)')
hold off
xlim([0 1000])

%% Write model file
ofile = sprintf('%s%s.model',odir,name);
fid = fopen(ofile,'w+');
fprintf(fid,'%s\n',name);
fprintf(fid,'%u %.4f %u\n',ifanis,tref,1);
fprintf(fid,'%u %u %u\n',N, Nk(1), Nk(2));
for kk = 1:N
    fprintf(fid,'%7.0f  %7.4f  %7.4f  %7.4f  %7.0f  %7.0f  %7.4f  %7.4f  %.4f \n',...
        R(kk),D(kk),Vpv(kk),Vsv(kk),Qk(kk),Qm(kk),Vph(kk),Vsh(kk),Eta(kk));
end
fclose(fid);