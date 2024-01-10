goodbye
% format long
load PWP_profile_3day
load PWP_met_ERA5
load timelim
%% plot the figures
figure('Position',[0 0 771 721])
figure('Position',[0 0 771 721])

% OutputFile='PWP_output_15Fluxes';
OutputFile='PWP_output_ControlFluxes_adv_fulldepth_2e_5';
load(OutputFile)

pres=gsw_p_from_z(-pwp_output.z,nanmean(profile.lat));
sa=gsw_SA_from_SP(pwp_output.s,pres,mean(profile.long),mean(profile.lat));
ct=gsw_CT_from_t(sa,pwp_output.t,pres);
ct_f=ct-gsw_CT_freezing(sa,pres);
rho0=gsw_rho(sa,ct,0)-1000;
[N2, p_mid]=gsw_Nsquared(sa,ct,pres,median(profile.lat));

day_start=met.day(1);
days=pwp_output.time+day_start;

% %% wind speed plot
% figure
% UU=met.U;
% Yixi_time_series(UU, met.day, 'Year Day', 'Wind Speed (m/s)')
% xlim([51.8271 105.7820]) % time(49:459)
% 
% % xlabel('')
% % set(gca,'XTickLabel',[])
% 
% UU_daily_movmean=movmean(UU,8,'omitnan');
% plot(met.day,UU_daily_movmean,'linewidth',1.7)

%% time series T
figure(1)
% h=subplot(2,1,1);
% set(h,'position',[0.1100 0.515 0.7750 0.41])
h=subplot(3,1,1);
set(h,'position',[0.1100 0.67 0.7750 0.2712])

plot_output_setup(OutputFile,'PWP_met_ERA5',ct_f,'\Theta above Freezing(ºC)')

caxis([-0.05 2.5]) 
ylim([0 600])
xlim([time_bgn time_end])

[C,h]=contour(days,pwp_output.z,rho0,[27.4:0.1:27.5],'LineWidth',2,'LineColor',[.3 .3 .3],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.3 .3 .3],'FontWeight','bold')
[C,h]=contour(days,pwp_output.z,rho0,[27 27.4],'LineWidth',2,'LineColor',[.79 .79 .79],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.79 .79 .79],'FontWeight','bold')

xlabel('')
set(gca,'XTick',55:5:105); set(gca,'XTickLabel',[])

% ax = gca;
% ax.XColor = 'w';
% ax.YColor = 'w'; 
% ax.GridColor = 'k';
% ax.MinorGridColor = 'k';
% set(gcf,'color','w')
% set(gca,'color','w')

%% time series S

figure(2)
h=subplot(3,1,1);
set(h,'position',[0.1100 0.67 0.7750 0.2712])

plot_output_setup(OutputFile,'PWP_met_ERA5',sa,'S_A (g kg^{-1})')
caxis([33.5 34.55])
xlim([time_bgn time_end])
ylim([0 600])
set(gca,'XTick',55:5:105); set(gca,'XTickLabel',[])
 xlabel('')

[C,h]=contour(days,pwp_output.z,rho0,[27 27.4:0.1:27.5],'LineWidth',2,'LineColor',[.3 .3 .3],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.3 .3 .3],'FontWeight','bold')
[C,h]=contour(days,pwp_output.z,rho0,[27 27.4],'LineWidth',2,'LineColor','k','ShowText','off');
clabel(C,h,'FontSize',17,'Color','k','FontWeight','bold')

% ax = gca;
% ax.XColor = 'w';
% ax.YColor = 'w'; 
% ax.GridColor = 'k';
% ax.MinorGridColor = 'k';
% set(gcf,'color','w')
% set(gca,'color','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for another model run
OutputFile='PWP_output_ControlFluxes';
load(OutputFile)
pres=gsw_p_from_z(-pwp_output.z,nanmean(profile.lat));
sa=gsw_SA_from_SP(pwp_output.s,pres,mean(profile.long),mean(profile.lat));
ct=gsw_CT_from_t(sa,pwp_output.t,pres);
ct_f=ct-gsw_CT_freezing(sa,pres);
rho0=gsw_rho(sa,ct,0)-1000;
[N2, p_mid]=gsw_Nsquared(sa,ct,pres,median(profile.lat));

day_start=met.day(1);
days=pwp_output.time+day_start;

%% time series T
figure(1)
h=subplot(3,1,2);
set(h,'position',[0.1100 0.37 0.7750 0.2712])

plot_output_setup(OutputFile,'PWP_met_ERA5',ct_f,'\Theta above Freezing(ºC)')

caxis([-0.05 2.5]) 
ylim([0 600])
xlim([time_bgn time_end])

[C,h]=contour(days,pwp_output.z,rho0,[27.5 27.5],'LineWidth',2,'LineColor',[.3 .3 .3],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.3 .3 .3],'FontWeight','bold')
[C,h]=contour(days,pwp_output.z,rho0,[27 27.4],'LineWidth',2,'LineColor',[.79 .79 .79],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.79 .79 .79],'FontWeight','bold')

xlabel('')
set(gca,'XTick',55:5:105); set(gca,'XTickLabel',[])

%% time series S
figure(2)
h=subplot(3,1,2);
set(h,'position',[0.1100 0.37 0.7750 0.2712])

plot_output_setup(OutputFile,'PWP_met_ERA5',sa,'S_A (g kg^{-1})')
caxis([33.5 34.55])
xlim([time_bgn time_end])
ylim([0 600])
set(gca,'XTick',55:5:105); set(gca,'XTickLabel',[])
xlabel('')

[C,h]=contour(days,pwp_output.z,rho0,[27 27.4:0.1:27.5],'LineWidth',2,'LineColor',[.3 .3 .3],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.3 .3 .3],'FontWeight','bold')
[C,h]=contour(days,pwp_output.z,rho0,[27 27],'LineWidth',2,'LineColor','k','ShowText','off');
clabel(C,h,'FontSize',17,'Color','k','FontWeight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for another model run
OutputFile='PWP_output_NoFlux';
load(OutputFile)
pres=gsw_p_from_z(-pwp_output.z,nanmean(profile.lat));
sa=gsw_SA_from_SP(pwp_output.s,pres,mean(profile.long),mean(profile.lat));
ct=gsw_CT_from_t(sa,pwp_output.t,pres);
ct_f=ct-gsw_CT_freezing(sa,pres);
rho0=gsw_rho(sa,ct,0)-1000;
[N2, p_mid]=gsw_Nsquared(sa,ct,pres,median(profile.lat));

day_start=met.day(1);
days=pwp_output.time+day_start;

%% time series T
figure(1)
h=subplot(3,1,3);
set(h,'position',[0.1100 0.07 0.7750 0.2712])

plot_output_setup(OutputFile,'PWP_met_ERA5',ct_f,'\Theta above Freezing(ºC)')

caxis([-0.05 2.5]) 
ylim([0 600])
xlim([time_bgn time_end])

[C,h]=contour(days,pwp_output.z,rho0,[27.4:0.1:27.5],'LineWidth',2,'LineColor',[.3 .3 .3],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.3 .3 .3],'FontWeight','bold')
[C,h]=contour(days,pwp_output.z,rho0,[27 27],'LineWidth',2,'LineColor',[.79 .79 .79],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.79 .79 .79],'FontWeight','bold')

xlabel('')
set(gca,'XTick',55:5:105); %set(gca,'XTickLabel',[])

%% time series S
figure(2)
h=subplot(3,1,3);
set(h,'position',[0.1100 0.07 0.7750 0.2712])

plot_output_setup(OutputFile,'PWP_met_ERA5',sa,'S_A (g kg^{-1})')
caxis([33.5 34.55])
xlim([time_bgn time_end])
ylim([0 600])
set(gca,'XTick',55:5:105); %set(gca,'XTickLabel',[])
xlabel('')

[C,h]=contour(days,pwp_output.z,rho0,[27 27.4:0.1:27.5],'LineWidth',2,'LineColor',[.3 .3 .3],'ShowText','off');
clabel(C,h,'FontSize',17,'Color',[.3 .3 .3],'FontWeight','bold')
[C,h]=contour(days,pwp_output.z,rho0,[27 27],'LineWidth',2,'LineColor','k','ShowText','off');
clabel(C,h,'FontSize',17,'Color','k','FontWeight','bold')




%% time series N2
% h=subplot(3,1,3);
% set(h,'position',[0.1100 0.07 0.7750 0.2712])
% 
% plot_output_setup(OutputFile,'PWP_met_ERA5',[nan(1,size(N2,2));N2],'log_1_0N^2 (rad^2 s^-^2)')
% caxis([-6 -3.6])
% ylim([0 600])
% xlim([time_bgn time_end])
% cb=colorbar('Ticks',-6:1:-4,'Position',[.89 .07 .015 .2712],'Color','w','FontSize',17);
% ylabel(cb,'log_1_0N^2 (rad^2 s^-^2)','FontSize',19,'color','w')
% 
% ax = gca;
% ax.XColor = 'w';
% ax.YColor = 'w'; 
% ax.GridColor = 'k';
% ax.MinorGridColor = 'k';
% set(gcf,'color','k')

% set(gcf, 'InvertHardcopy', 'off')
