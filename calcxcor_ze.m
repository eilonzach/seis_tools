% calcxcor.m   calculate cross-correlation corrections to arrival times
% This script uses assoc, arrival, and predarr tables (as well as site and orid)
%
% The predarr table contains the predicted arrival times and slownesses. If
% this does not exist it must be made using dbassoc2predarr.m, which
% requires assoc and arrival tables (as well as origin and site)
close all %, clear all
% STARTS - SET THESE
dbnam = 'pngbodydb';
dbdir = '/Users/zeilon/Documents/MATLAB/PNG_tomog/PNGbodydb/';
odir = '/Users/zeilon/Documents/MATLAB/PNG_tomog/';

iphase='S';
chan='BHN';
%times (s) for extraction window
pretime = 75.;     
posttime = 75.;
%time window for cross-correlation
% prex=2.;   
% postx=5.;   
taperx = 0.1;
% DEFAULT hicut, locut filter corners for xcor - [n x 1] column vectors of n filters
fhi = [1]';      
flo = [0.05]';
npoles = [4]';
% max allowed xcor-lag to test; CHECK
lagmax=8; 

imode=1;      % aligns by picked arrival

startorid = 1; 

%% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX %%
%% XXXXXXXXXXXXXXXXXXXXXXXXX  REST DOES WORK  XXXXXXXXXXXXXXXXXXXXXXXXXX %%
%% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX %%
ofile = strcat(odir,'xcor_results_',iphase,'_',chan(3));
ofile = strcat(odir,'xcor_resultsJUNK_',iphase,'_',chan(3));

cd(dbdir);
db=dbopen(dbnam,'r+');
dbar=dblookup_table(db,'arrival');
dbas=dblookup_table(db,'assoc');
dbor=dblookup_table(db,'origin');
dbsi=dblookup_table(db,'site');
dbpa=dblookup_table(db,'predarr');

orids=unique(dbgetv(dbas,'orid'));
norids=length(orids);

allstas = dbgetv(dbsi,'sta');

figure(1), set(gcf,'Position',[600   50 1300 930])
figure(2), set(gcf,'Position',[2000 250 1320 650])

%% To do only a subset of orids and put in its own result structure
% orids = nprobors;
% norids=length(orids);

%% Results structures
% results in [nsta x norid] blocks, where norid is all the successful (and
% hence saved) xcorrelation measurements. NaN values for stations with no
% measurements for a given orid
ovrw = 1; % default is make new results structures
if exist(strcat(ofile,'.mat'),'file')==2
ovrw=menu('WARNING, some results already saved... ','Overwrite and restart','Append to existing');
end
if ovrw==1;
alltabs   = NaN(size(allstas,1),0);
alldcor   = NaN(size(allstas,1),0);
alldcstd   = NaN(size(allstas,1),0);
allacor   = NaN(size(allstas,1),0);
allseazs  = NaN(size(allstas,1),0);
allgcarcs = NaN(size(allstas,1),0);
allslows  = NaN(size(allstas,1),0);
allarids  = NaN(size(allstas,1),0);
alltrefs  = NaN(size(allstas,1),0);
alltprds  = NaN(size(allstas,1),0);
allxcor_parms = struct('orid',0,'evdata',struct([]),'chan',chan,'filter',[0 0 0],'taperx',taperx,'window',[0 0],'lagmax',[lagmax],'prepost',[pretime posttime],'note','');
xcor_results = struct([]);
elseif ovrw==2
    load(strcat(ofile,'.mat'));
end

%% Get going...
while ~any(orids==startorid)
    startorid = startorid+1;
end

for ie = find(orids==startorid):norids
orid = orids(ie)
figure(1); clf
try
  [data,traces,tts,nsamps,samprates,stas,chans,gcarcs,seazs,esazs,trefs,arids,evdata]...
 = dbgetdata(db,'evt',orid,'p', iphase,'c', chan,'w',[pretime posttime],'mode',imode);
catch
    continue
end

% get predicted times and slownesses
[xarids,ttpreds,slows] = dbgetv(dbpa,'arid','time','slow');
ttpreds = ttpreds(ismember(xarids,arids));
slows = slows(ismember(xarids,arids)); 
% put predicted times into station order
[~,junk]=sort(arids); [~,sta_arid_order] = sort(junk);
ttpreds = ttpreds(sta_arid_order);
slows = slows(sta_arid_order);

%% get rest of info about output
nsta = length(stas);
npt = size(traces,1);
nchan = length(chans);
samprate = unique(samprates);
dt = 1./samprate;

%% START PROCESSING THIS EVENT
nextorid = 0;
delsta = '';
while nextorid==0;
%% If a station is to be deleted...
if ~isempty(delsta)
delind = listdlg('Liststring',stas); 
data(delind,:) = [];
traces(:,nchan*(delind-1)+(1:nchan))=[];
tts(:,nchan*(delind-1)+(1:nchan))=[];
nsamps(delind)=[];  samprates(delind)=[];   stas(delind)=[];
gcarcs(delind)=[];  seazs(delind)=[];       esazs(delind)=[]; 
trefs(delind)=[];   arids(delind)=[];       slows(delind)=[];
ttpreds(delind)=[];
nsta = length(stas);
npt = size(traces,1);
end

%% plot for windowing
figure(1); clf
for is=1:nsta
hp1 = subplot(131); hold on
plot((0:npt-1)*dt - pretime,traces(:,is)./max(abs(traces(:,is)))+is,'Linewidth',1.5);
ylim([0 is+1]); xlim([-pretime posttime-dt])
set(gca,'YTick',1:nsta,'YTickLabel',stas,'FontWeight','bold')
title(sprintf('Original data\nORID %u, seaz %.0f, gcarc %.0f\ndep %.0fkm, Ms %.1f',orid,mean(seazs),mean(gcarcs),evdata.edep,evdata.mag),'Fontsize',18)
hold off
end

%% Option to edit filter
flt = inputdlg({'Enter  max  period','Enter  min  period','Enter  N  filter  poles'},...
                   'Filter details',1,{num2str(1/flo),num2str(1/fhi),num2str(npoles)});
flt = str2num(char(flt));
fhi = 1/flt(2);
flo = 1/flt(1);
npoles = flt(3);

cont=0;
while cont==0
try delete(hp2); end
figure(1); subplot(131)
[wx,~] = ginput(2); 
wx = sort(round_level(wx,2*dt));
try delete(h1); end
hold on; h1 = plot([1;1]*wx',[0;nsta+1]*[1,1],'--k','Linewidth',1.5); hold off
if wx(2)-wx(1)<1/flo
uiwait(warndlg('Window  length  shorter  than longest  period  in  filter... curtailing longest period','Window nüdge','modal'))
end

%% Treat data
cleaning_parm = struct('samprate',samprate,'pretime',pretime,'prex',-wx(1),'postx',wx(2),...
                       'taperx',taperx,'fhi',fhi,'flo',flo,'npoles',npoles,'norm',1);
[ datwf,datf,datc,fdeets,ttsf,ttsc ] = data_clean( traces,cleaning_parm );
fhi = fdeets(1);
flo = fdeets(2);

if any(isnan(datwf))
uiwait(warndlg('Broken filter... try fewer poles','Filter fail','modal'))
cont=1;
continue
end

%% Plot windowed and treated data
figure(1)
for is=1:nsta
hp2 = subplot(132); hold on
plot(ttsc, datf(:,is)./max(abs(datf(:,is)))+is,'Linewidth',1);
plot(ttsf, datwf(:,is)./max(abs(datwf(:,is)))+is,'r','Linewidth',1.5);
ylim([0 is+1]); xlim([ttsc(1) ttsc(end)])
set(gca,'YTick',[])
title('Cleaned, filtered, windowed','Fontsize',18)
hold off
end % loop on stas

%% Decide to re-window
uic1 = uicontrol('Style', 'pushbutton', 'String', 'OK',...
        'Position', [1000 610 85 20],'Callback', 'cont=1; uiresume(1)');
uic2 = uicontrol('Style', 'pushbutton', 'String', 'Re-window',...
        'Position', [1100 610 85 20],'Callback', 'cont=0;uiresume(1)');
uiwait(1)
% waitforbuttonpress; pause(0.1)
delete(uic1); delete(uic2);
end   

% go back to filter choice if filter is broken
if any(isnan(datwf))
continue
end

%% Do xcorr
figure(2); 
disp('CROSS-CORRELATING...');
[dcor, dcstd,~,acor]=xcortimes(datwf, dt, -ttsf(1),lagmax,1);  % Subtracted from timebase
%  will add dcor to the "actual" picks
  title('cross-correlated records')
  xlabel('lag time, s')
  ht = labelText(gcf); set(ht,'FontSize',18);
%   nscc=length(ttsf);
%   stak=zeros(nscc,1);
%   tstak=reshape(ttsf,nscc,1);
%   for is=1:nsta
%     stak=stak+interp1(tstak-dcor(is),datwf(:,is),tstak);
%   end
%   stak=stak./nsta;
%   jmx=nscc-ceil(max(dcor)./dt);
%   jmn=ceil(abs(min(dcor))./dt)+1;
%   ida=find(isfinite(stak) & (1:nscc)'<jmx & (1:nscc)'>jmn );
%   plot(tstak,stak,'r','LineWidth',2);
%   clear('acor')
%   
%   for is=1:nsta
%     acor(is)=sum(stak(ida).*datwf(ida+round(dcor(is)./dt),is))./(std(stak(ida),1).*std(datwf(ida+round(dcor(is)./dt),is),1))./length(ida);
%   end
  

%% plot aligned traces with more tail
cleaning_parm.prex = pretime/2;%min([-2*wx(1),pretime]);
cleaning_parm.postx = posttime/2;%min([2*wx(2),posttime-dt]);
[ datwf,~,~,~,ttsf,~ ] = data_clean( traces,cleaning_parm );

figure(1)
hp3 = subplot(133); cla; hold on
for is=1:nsta
plot(ttsf-dcor(is),datwf(:,is)./max(abs(datwf(:,is)))+is,'Linewidth',1.5);
ht = text(ttsf(end)+0.05*(ttsf(end)-ttsf(1)),is,sprintf('acor = %.3f\n  dcor = %.3f',acor(is),dcor(is)));
set(ht,'FontSize',13,'FontWeight','bold');
end
hold on; plot([0;0],[0;nsta+1]*[1,1],'--k','Linewidth',1); hold off
ylim([0 is+1]); xlim([ttsf(1) ttsf(end)])
set(gca,'YTick',[])
title(sprintf('Mean acor = %.3f\nAligned & cleaned etc.',mean(acor)),'Fontsize',18)
hold off

%% Decide to re-window
delsta = '';
uic1 = uicontrol('Style', 'pushbutton', 'String', 'Save and move on',...
        'Position', [350 30 120 30],'Callback', 'nextorid=1;');
uic2 = uicontrol('Style', 'pushbutton', 'String', 'Discard and move on',...
        'Position', [500 30 120 30],'Callback', 'nextorid=2;');
uic3 = uicontrol('Style', 'pushbutton', 'String', 'Re-do',...
        'Position', [650 30 120 30],'Callback', 'nextorid=0;');
uic4 = uicontrol('Style', 'pushbutton', 'String', 'Re-do and ignore a sta',...
        'Position', [800 30 140 30],'Callback', 'nextorid=0; delsta = ''y'';');
waitforbuttonpress; pause(0.2)
delete(uic1); delete(uic2); delete(uic3); delete(uic4);
end % while doing each station loop
return
%% SAVE
if nextorid==1
alltabs(:,size(alltabs,2)+1) = NaN;
alltabs(ismember(allstas,stas),size(alltabs,2)) = tts(1,:)'+pretime+dcor; 
alldcor(:,size(alltabs,2)) = NaN;
alldcor(ismember(allstas,stas),size(alltabs,2)) = dcor;                     
alldcstd(:,size(alltabs,2)) = NaN;
alldcstd(ismember(allstas,stas),size(alltabs,2)) = dcstd;                     
allacor(:,size(alltabs,2)) = NaN;
allacor(ismember(allstas,stas),size(alltabs,2)) = acor;                     
allseazs(:,size(alltabs,2)) = NaN;
allseazs(ismember(allstas,stas),size(alltabs,2)) = seazs;                   
allgcarcs(:,size(alltabs,2)) = NaN;
allgcarcs(ismember(allstas,stas),size(alltabs,2)) = gcarcs;                 
allslows(:,size(alltabs,2)) = NaN;
allslows(ismember(allstas,stas),size(alltabs,2)) = slows; 
allarids(:,size(alltabs,2)) = NaN;
allarids(ismember(allstas,stas),size(alltabs,2)) = arids;                   
alltrefs(:,size(alltabs,2)) = NaN;
alltrefs(ismember(allstas,stas),size(alltabs,2)) = trefs;
alltprds(:,size(alltabs,2)) = NaN;
alltprds(ismember(allstas,stas),size(alltabs,2)) = ttpreds;
xcor_parms = struct('orid',orid,'evdata',evdata,'chan',chans,'filter',[fhi flo npoles],'taperx',taperx,'window',wx','lagmax',lagmax,'prepost',[pretime posttime],'note','');
allxcor_parms(1,size(alltabs,2)) = xcor_parms;

xcor_results = struct('alltabs',alltabs,'alldcor',alldcor,'alldcstd',alldcstd,...
                      'allacor',allacor,'allseazs',allseazs,...
                      'allgcarcs',allgcarcs,'allslows',allslows,...
                      'allarids',allarids,'alltrefs',alltrefs,...
                      'alltprds',alltprds,'allxcor_parms',allxcor_parms);
save(ofile,'-struct','xcor_results','alltabs','alldcor','alldcstd',...
                      'allacor','allseazs','allgcarcs','allslows',...
                      'allarids','alltrefs','alltprds','allxcor_parms');
save(ofile,'allstas','-append');
end

%% Check alignment on feature
% figure(3); clf
% for is = 1:nsta
% try
%   [data,trace2,tt2] = dbgetdata(db,'time',alltabs(end,strcmp(allstas,stas(is)))+[-pretime posttime],'c',chan,'s',char(stas(is)));
% catch; continue
% end
% hold on
% plot(dt*(0:(npt-1)) - pretime,trace2/max(abs(trace2))+is)
% hold off
% end
% set(gca,'YTick',1:nsta,'YTickLabel',flipud(stas),'FontWeight','bold')


end % loop on orids
dbclose(db);
cd(odir);
  
  
  