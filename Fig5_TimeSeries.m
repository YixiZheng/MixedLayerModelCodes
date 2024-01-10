bye
load EM959_calib

spelim=long<-111 & lat>-74.2;
time_bgn=min(time(spelim));
time_end=max(time(spelim));

ct=ct(spelim);sa=sa(spelim);
time=time(spelim);dep=dep(spelim);pres=pres(spelim);
long=long(spelim);lat=lat(spelim);
clearvars spelim 

%% daily rho
load EM959_binned_calib_2m
rho0_b=gsw_rho(sa_b,ct_b,0)-1000;
j=0;
time_daily=floor(min(time_b)):ceil(max(time_b))-1;
for i=time_daily
    j=j+1;
    rho0_daily(j,:)=nanmean(rho0_b(time_b>=i & time_b<i+1,:));
end

%% time series T
figure('position',[759 176 1170 1004])

h=subplot(3,1,1);
set(h,'position',[0.1100 0.67 0.7750 0.2712])

ct_f=ct-gsw_CT_freezing(sa,pres);
Yixi_time_series_scatter(time,-dep,13,ct_f,'Year','\Theta above Freezing(ºC)')
ylabel('Depth (m)','fontsize',23)
caxis([-0.05 2.5]) 
% caxis([0 4]) 
xlim([min(time) max(time)])
ylim([0 600])
%ylim([0 1100])

load obs_MLD % from obs_FindMLD.m

plot(obs_MLD_time,obs_MLD,'LineWidth',3,'color','k')
plot(obs_MLD_time,obs_MLD,'LineWidth',1,'color','w')

[C,h]=contour(time_daily,-dep_b,rho0_daily',[27 27.4:0.1:27.5],'LineWidth',2,'LineColor',[.3 .3 .3],'ShowText','on');
clabel(C,h,'FontSize',17,'Color',[.3 .3 .3],'FontWeight','bold')
[C,h]=contour(time_daily,-dep_b,rho0_daily',[27 27],'LineWidth',2,'LineColor',[.79 .79 .79],'ShowText','on');
clabel(C,h,'FontSize',17,'Color',[.79 .79 .79],'FontWeight','bold')
set(gca,'XTick',55:5:105); set(gca,'XTickLabel',[])
xlabel('')
cb=colorbar('Ticks',0:1:2.5,'Position',[.89 .67 .015 .2712],'Color','k','FontSize',23);
ylabel(cb,'\Theta above Freezing(ºC)','FontSize',23,'color','k')


% ax = gca;
% ax.XColor = 'w';
% ax.YColor = 'w'; 
% ax.GridColor = 'k';
% ax.MinorGridColor = 'k';
% set(gcf,'color','k')

%% time series S
h=subplot(3,1,2);
set(h,'position',[0.1100 0.37 0.7750 0.2712])

Yixi_time_series_scatter(time,-dep,13,sa,'Year Day','S_A (g kg^-^3)')
ylabel('Depth (m)','fontsize',23)
% caxis([34 35])
caxis([33.5 34.55])
xlim([min(time) max(time)])
%ylim([0 1100])
ylim([0 600])

load obs_MLD

plot(obs_MLD_time,obs_MLD,'LineWidth',3,'color','k')
plot(obs_MLD_time,obs_MLD,'LineWidth',1,'color','w')

[C,h]=contour(time_daily,-dep_b,rho0_daily',[27 27.4:0.1:27.5],'LineWidth',2,'LineColor',[.3 .3 .3],'ShowText','on');
clabel(C,h,'FontSize',17,'Color',[.3 .3 .3],'FontWeight','bold')
% [C,h]=contour(time_daily,-dep_b,rho0_daily',[27 27],'LineWidth',1,'LineColor',[.99 .99 .99],'ShowText','on');
% clabel(C,h,'FontSize',17,'Color',[.99 .99 .99])
set(gca,'XTick',55:5:105); set(gca,'XTickLabel',[])
xlabel('')

cb=colorbar('Ticks',33.5:0.5:34.5,'Position',[.89 .37 .015 .2712],'Color','k','FontSize',23);
ylabel(cb,'S_A (g kg^-^3)','FontSize',23,'color','k')



% ax = gca;
% ax.XColor = 'w';
% ax.YColor = 'w'; 
% ax.GridColor = 'k';
% ax.MinorGridColor = 'k';
% set(gcf,'color','k')

%% time series N2
h=subplot(3,1,3);
set(h,'position',[0.1100 0.07 0.7750 0.2712])

[N2, p_mid]=gsw_Nsquared(sa,ct,pres,median(lat));
N2(abs(diff(pres))>10)=NaN;
N2(abs(N2)==Inf)=NaN;
Yixi_time_series_scatter(time(1:end-1)+diff(time)/2,-(dep(1:end-1)+diff(dep)/2),13,log10(N2),'Year Day','log_1_0N^2 (rad^2 s^-^2)')
ylabel('Depth (m)','fontsize',23)
caxis([-6 -3.6])

load obs_MLD

plot(obs_MLD_time,obs_MLD,'LineWidth',3,'color','k')
plot(obs_MLD_time,obs_MLD,'LineWidth',1,'color','w')

[C,h]=contour(time_daily,-dep_b,rho0_daily',[27 27.4:0.1:27.5],'LineWidth',2,'LineColor','k','ShowText','on');
clabel(C,h,'FontSize',17,'Color','k','FontWeight','bold')

cb=colorbar('Ticks',-6:1:-4,'Position',[.89 .07 .015 .2712],'Color','k','FontSize',23);
ylabel(cb,'log_1_0N^2 (s^-^2)','FontSize',23,'color','k')

xlim([min(time) max(time)])
%ylim([0 1100])
ylim([0 600])


% ax = gca;
% ax.XColor = 'w';
% ax.YColor = 'w'; 
% ax.GridColor = 'k';
% ax.MinorGridColor = 'k';
% set(gcf,'color','k')
% 
% set(gcf, 'InvertHardcopy', 'off')