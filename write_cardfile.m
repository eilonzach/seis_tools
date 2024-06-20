function card = write_cardfile( ofile,Z,vpv,vsv,rho,Qk,Qmu,vph,vsh,eta,tref,minlaydz ) %#ok<INUSL>
% card = write_cardfile( ofile,Z,vp,vs,rho,[Qk=prem_value],[Qmu=prem_value],[vph=vpv],[vsh=vsv],[eta=1],[tref=1],[minlaydz=20])
%   
% This function writes a cardfile - format to be used for (e.g.) MINEOS.
% The function expects that you provide vp,vs,rho at a set of upper mantle
% depths. Below the deepest velocity provided, the function will smoothly
% transition into anisotropic PREM


if isempty(ofile)
    ifwrite=false;
else
    ifwrite=true;
end
if nargin < 8 || isempty(vph) 
    vph = vpv;
end
if nargin < 9 || isempty(vsh) 
    vsh = vsv;
end
if nargin < 10 || isempty(eta) 
    eta = 1*ones(size(Z));
end
if nargin < 11 || isempty(tref) 
    tref = 1;
end
if nargin < 12 || isempty(minlaydz) 
    minlaydz = 20;
end

% if tref > 0; 
%     tref = - tref; 
%     fprintf('Making mineos attenuation reference velocity negative for Russell Mineos code. '); 
%     % brb20240606 TODO add option to allow positive if using Zach's version of Mineos. 
% end

global prem_anisotropic prem_isotropic
if isempty(prem_isotropic)
    prem_isotropic = prem;
end
if isempty(prem_anisotropic)
    prem_anisotropic = prem_perfect('SPVW',0.25);
end
    
Re = 6371;
gradz = 100; % thickness, in km, of region of linear grading between my model and PREM
% minlaydz = 20;

%% Fields
flds = {'rho','vpv','vsv','vph','vsh','Qk','Qmu','eta'};

%% Grab PREM
if isequal(vpv,vph) && isequal(vsv,vsh) && all(eta==1) % if totally isotropic
    prem_mod = prem_isotropic;
    % refine depths
    zz = prem_mod.depth;
    while any(diff(zz)>minlaydz)
        irep = find(diff(zz)>minlaydz);
        zz = sort([zz;mean([zz(irep),zz(irep+1)],2)]);
    end
    
    prem_mod = prem('depths',zz,'dcbelow',false);
    prem_mod.vpv = prem_mod.vp; prem_mod.vph = prem_mod.vp;
    prem_mod.vsv = prem_mod.vs; prem_mod.vsh = prem_mod.vs;
    prem_mod.eta = ones(size(prem_mod.depth));
else % anisotropic!
    prem_mod = prem_anisotropic; % if anisotropic - slower!
end
prem_mod.Qk = prem_mod.qk;
prem_mod.Qmu = prem_mod.qu;

%% Q values
if nargin<6 || isempty(Qk)
    Qk = linterp(prem_mod.depth,prem_mod.Qk,Z); %#ok<NASGU>
end
if nargin<7 || isempty(Qmu)
    Qmu = linterp(prem_mod.depth,prem_mod.Qmu,Z); %#ok<NASGU>
end

%% Work out bottom of model + stitching details
maxz = max(Z);
igrad = (prem_mod.depth>maxz) & (prem_mod.depth<= maxz+gradz); % indices of depths in gradation region
fgrad = (maxz+gradz-prem_mod.depth(igrad))/gradz; %#ok<NASGU> % gradation factor (between 1 at the top and 0 at the bottom of gradation region)
bgrad = prem_mod.depth> maxz+gradz;

%% Make the collated vectors for each datatype
for ii = 1:length(flds)
    eval(sprintf('dval = %s(end) - linterp(prem_mod.depth,prem_mod.%s,maxz);',flds{ii},flds{ii}))
    eval(sprintf('card.%s = [%s;dval*fgrad + prem_mod.%s(igrad);prem_mod.%s(bgrad)];',flds{ii},flds{ii},flds{ii},flds{ii}))
end
card.depth = [Z;prem_mod.depth(igrad);prem_mod.depth(bgrad)];
card.R = Re - card.depth;
card.Qmu(isinf(card.Qmu) & card.depth>2500) = 0 ; % fix issue where linterp makes zero Qmu inf.
card.Qmu(isnan(card.Qmu) & card.depth<10) = 999 ; % fix issue where linterp makes Inf Qmu nan.
card.Qmu(isinf(card.Qmu)) = 999 ; % fix issue where linterp makes Inf Qmu nan.

%% make sure in units of km (or km/s, or kg/m^2)
if max(card.R)>7000, card.R = card.R/1000; end
if max(card.rho)>100, card.rho = card.rho/1000; end
if max(card.vpv)>100, card.vpv = card.vpv/1000; end
if max(card.vph)>100, card.vph = card.vph/1000; end
if max(card.vsv)>100, card.vsv = card.vsv/1000; end
if max(card.vsh)>100, card.vsh = card.vsh/1000; end

%% find node numbers
N = length(card.R);
ocind = find(card.Qmu==0 & card.vsv==0);
Nic = N-ocind(end); % number of inner core nodes
Noc = length(ocind)+Nic;

%% Write model file
if ifwrite
% edit cardfile name so no periods
ofile_print = ofile;
ofile_print(regexp(ofile_print,'\.')) = '_';

fid = fopen(ofile,'w+');
fprintf(fid,'%s\n',ofile_print);
fprintf(fid,'  %u   %.4f   %u\n',1,tref,1); % make anis flag = 1 even if not truly anisotropic!
fprintf(fid,'  %u   %u   %u\n',N, Nic, Noc);
for ii = 1:N
    kk = N+1 - ii; % reverse order to have zero radius first
    fprintf(fid,'%7.0f.%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n',...
        1000*card.R(kk),1000*card.rho(kk),1000*card.vpv(kk),1000*card.vsv(kk),card.Qk(kk),card.Qmu(kk),1000*card.vph(kk),1000*card.vsh(kk),card.eta(kk));
end
fclose(fid);
fprintf('cardfile %s written\n',ofile);
end

end

