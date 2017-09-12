function plot_KERNELS( phV_kernels,ifanis,ax1,ax2,ofile )
%plot_KERNELS( phV_kernels,ifanis=false,ax1=[],ax2=[] )


if nargin<2 || isempty(ifanis)
    ifanis = false;
end
    
if nargin < 3 || isempty(ax1) || isempty(ax2)
    figure(88), clf; set(gcf,'pos',[680 385 999 713]);
    ax1 = subplot(1,2,1); cla, hold on;
    ax2 = subplot(1,2,2); cla, hold on;   
end

if nargin < 5|| isempty(ofile), ifsave = false; else ifsave = true; end

np = length(phV_kernels);
swperiods = zeros(np,1); for ip = 1:np, swperiods(ip) = phV_kernels{ip}.period; end



%% plot kernels
h = zeros(np,1);
for ip = 1:np
    if ifanis
    h(ip)=plot(ax1,phV_kernels{ip}.Vsv,phV_kernels{ip}.Z/1e3,'linewidth',2);
        plot(ax1,phV_kernels{ip}.Vsh,phV_kernels{ip}.Z/1e3,':','linewidth',2)
        plot(ax2,phV_kernels{ip}.Vpv,phV_kernels{ip}.Z/1e3,'linewidth',2)
        plot(ax2,phV_kernels{ip}.Vph,phV_kernels{ip}.Z/1e3,':','linewidth',2)
        labstrs = 'SV=solid, SH=dash'; labstrp = 'PV=solid, PH=dash';
    else
    h(ip)=plot(ax1,phV_kernels{ip}.Vsv+phV_kernels{ip}.Vsh,phV_kernels{ip}.Z/1e3,'linewidth',2);
        plot(ax2,phV_kernels{ip}.Vpv+phV_kernels{ip}.Vph,phV_kernels{ip}.Z/1e3,'linewidth',2)        
        labstrs = ''; labstrp = '';
    end
end
%% axes properties
set(ax1,'ydir','reverse','ylim',[0 400],'fontsize',16)
set(ax2,'ydir','reverse','ylim',[0 400],'fontsize',16)
xlabel(ax1,'Sensitivity to $V_S$',...
    'interpreter','latex','fontsize',22)
title(ax1,labstrs);
xlabel(ax2,'Sensitivity to $V_P$',...
    'interpreter','latex','fontsize',22)
title(ax2,labstrp);
ylabel(ax1,'Depth (km)','interpreter','latex','fontsize',22)
ylabel(ax2,'Depth (km)','interpreter','latex','fontsize',22)
%% legend
hl = legend(ax1,h,num2str(round(swperiods)),'location','southeast');
set(get(hl,'title'),'string','Period (s)','fontsize',18,'fontweight','bold')
set(hl,'fontsize',16,'Linewidth',2)

if ifsave
    save2pdf(gcf,ofile);
end



end

