function [ hfig,ax ] = plot_cardfile( cardfile )
% [ hfig,ax ] = plot_cardfile( cardfile )

% read in
[ model ] = read_cardfile( cardfile );

%% set up plot
hfig = figure(55); clf, set(hfig,'pos',[17 54 1320 895]);

% different axes
ax1 = axes('pos',[0.07 0.08 0.32 0.88]); hold on
ax2 = axes('pos',[0.60 0.23 0.32 0.73]); hold on
ax3 = axes('pos',[0.42 0.08 0.15 0.88]); hold on

% plot velocities and density
plot(ax1,model.Vpv,model.Z,'-b',model.Vph,model.Z,'--b','linewidth',2)
plot(ax2,model.Vpv,model.Z,'-b',model.Vph,model.Z,'--b','linewidth',2)

plot(ax1,model.Vsv,model.Z,'-r',model.Vsh,model.Z,'--r','linewidth',2)
plot(ax2,model.Vsv,model.Z,'-r',model.Vsh,model.Z,'--r','linewidth',2)

plot(ax1,model.rho,model.Z,'-k','linewidth',2)
plot(ax2,model.rho,model.Z,'-k','linewidth',2)

plot(ax3,model.Qk,model.Z,'-k',model.Qm,model.Z,'-m','linewidth',2)


% pretty
set([ax1,ax2,ax3],'ylim',[0 max(model.Z)],'ydir','reverse',...
       'linewidth',2,'box','on','layer','top',...
       'fontsize',16)
set(ax2,'ylim',[0 700],'yaxislocation','right');
set(ax3,'xlim',[0 700],'yticklabel','');
   
ylabel(ax1,'Depth (km)','fontsize',20)
ylabel(ax2,'Depth (km)','fontsize',20)
xlabel(ax1,'Vp(b), Vs(r), rho(k)','fontsize',20)
xlabel(ax2,'Vp(b), Vs(r), rho(k)','fontsize',20)
xlabel(ax3,'Qmu','fontsize',20)


ax = [ax1,ax2,ax3];



end