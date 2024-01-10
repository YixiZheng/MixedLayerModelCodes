bye
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just to confirm, Yixi is still using this code to run PWP, 30/Oct/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just to confirm, Yixi is still using this code to run PWP, 17/May/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just to confirm, Yixi is still using this code to run PWP, 21/Dec/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this one is a running experiment for the 22/Jul forward manuscript. Yixi
% generate this version because the previous PWP_run.m is really messy. At
% the same time Yixi doesn't want to spend too much time editing it, while
% trying not to delete useful things. So Yixi just made a new version,
% based on the old one. Give her mercy, she's getting crazy because of her thesis

% this version does only one thing, running all experiment in one time.
% so far (22 Jul 2022), it makes sense
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this one is a totally new experiment, designed to test the influence of ice coverage to
% the mixed layer deepening.
%
% the Amelt = the area of ice, which can cause brine rejection and heat
% consumption when sea ice is melting, can understand as the "percentage of
% extra heat that will be used to form sea ice"
% the Agrow = the area of ice, which can cause brine rejection and heat
% consumption when sea ice is growing
% the A = the area of ice, which can block the heat exchange and momentum
% exchange
%
% the hypothesis that Yixi has is that, if the wind keeps blowing the sea
% ice away, the A and Amelt/Agrow are going to be different - i.e. the sea
% ice keeps growing, and rejecting brine and consuming heat, but at the
% same time, the A is nearly zero, which means that we still have strong
% heat loss from the ocean to the air.
%
% so far (27 Feb 2022), it makes sense
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% for testing different grids
% mPWP_phys_EM959('PWP_met_ERA5_1', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_grid1',1,1,0) 
% % mPWP_phys_EM959('PWP_met_ERA5_1', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot25_grid1',1,1,0.25) 
% % mPWP_phys_EM959('PWP_met_ERA5_1', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot5_grid1',1,1,0.5) 
% % mPWP_phys_EM959('PWP_met_ERA5_1', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot75_grid1',1,1,0.75) 
% % mPWP_phys_EM959('PWP_met_ERA5_1', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_grid1',1,1,1)
% % 
% mPWP_phys_EM959('PWP_met_ERA5_2', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_grid2',1,1,0) 
% % mPWP_phys_EM959('PWP_met_ERA5_2', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot25_grid2',1,1,0.25) 
% % mPWP_phys_EM959('PWP_met_ERA5_2', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot5_grid2',1,1,0.5) 
% % mPWP_phys_EM959('PWP_met_ERA5_2', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot75_grid2',1,1,0.75) 
% % mPWP_phys_EM959('PWP_met_ERA5_2', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_grid2',1,1,1)
% % 
% mPWP_phys_EM959('PWP_met_ERA5_3', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_grid3',1,1,0) 
% % mPWP_phys_EM959('PWP_met_ERA5_3', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot25_grid3',1,1,0.25) 
% % mPWP_phys_EM959('PWP_met_ERA5_3', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot5_grid3',1,1,0.5) 
% % mPWP_phys_EM959('PWP_met_ERA5_3', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot75_grid3',1,1,0.75) 
% % mPWP_phys_EM959('PWP_met_ERA5_3', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_grid3',1,1,1)
% 
% mPWP_phys_EM959('PWP_met_ERA5_4', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_grid4',1,1,0) 
% % mPWP_phys_EM959('PWP_met_ERA5_4', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot25_grid4',1,1,0.25) 
% % mPWP_phys_EM959('PWP_met_ERA5_4', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot5_grid4',1,1,0.5) 
% % mPWP_phys_EM959('PWP_met_ERA5_4', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot75_grid4',1,1,0.75) 
% % mPWP_phys_EM959('PWP_met_ERA5_4', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_grid4',1,1,1)

%% for testing heat flux 
% mPWP_phys_EM959('PWP_met_ERA5_hourly_north', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_north',1,1,0) 

%% for testing heat flux 
% note, when the Aicecover=0.25, the fluxes are decreased by 0.25, so the fluxes are in the 0.75 of their original value
% for the nagetive values, it means the fluxes are actaully increased by, say, 0.25, so will be positive
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes',1,1,0) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_075fluxes',1,1,0.25) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_05fluxes',1,1,0.5) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_025fluxes',1,1,0.75) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_NoFlux',1,1,1) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_125Fluxes',1,1,-.25) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_15Fluxes',1,1,-.5)
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_175Fluxes',1,1,-.75)
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_DoubleFluxes',1,1,-1)
%% for testing momentum flux
% % mPWP_phys_EM959_onlyMomentum('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_onlyMomentum',1,1,0) 
% % mPWP_phys_EM959('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_ControlMomentum',1,1,0) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_0tau_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_NoMomentum',1,1,0) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_0tau_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_NoFlux_NoMomentum',1,1,1) 

mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_0tau_0U_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_NoWind',1,1,0) 
mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5_hourly_0tau_0U_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_NoFlux_NoWind',1,1,1) 
% mPWP_phys_EM959('PWP_met_ERA5_025tau', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_025tau',1,1,0) 
% mPWP_phys_EM959('PWP_met_ERA5_05tau', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_05tau',1,1,0) 
% mPWP_phys_EM959('PWP_met_ERA5_075tau', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_075tau',1,1,0) 
% mPWP_phys_EM959('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_1tau',1,1,0) 
%% for testing freshwater flux
% mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_nofw',1,1,0) 
% mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot25_nofw',1,1,0.25) 
% mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot5_nofw',1,1,0.5) 
% mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_dot75_nofw',1,1,0.75) 
% mPWP_phys_EM959_OnlyHeatBlocked('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_1_nofw',1,1,1)

%% for advection
mPWP_phys_EM959_advection_onlyHeatBlocked('PWP_met_ERA5_hourly_albedoed', 'PWP_profile_3day_EOS', 'PWP_output_ControlFluxes_adv_fulldepth_2e_5',1,1,0,-2e-5,-2e-5)

%% plot the figures
% load PWP_profile_3day_EOS
% load PWP_met_ERA5_0tau
% OutputFile='PWP_output_1_1_0';
load(OutputFile)

pres=gsw_p_from_z(-pwp_output.z,nanmean(profile.lat));
sa=gsw_SA_from_SP(pwp_output.s,pres,mean(profile.long),mean(profile.lat));
ct=gsw_CT_from_t(sa,pwp_output.t,pres);
ct_f=ct-gsw_CT_freezing(sa,pres);

day_start=met.day(1);
time=pwp_output.time+day_start;
pres_plot=repmat(pres,1,length(time));

%% heat content plot
figure
hc=heat_content(sa',ct',pres_plot',mean(diff(pwp_output.z)),67);
% hc=heat_content(sa',ct',pres_plot',mean(diff(pwp_output.z)),200);
Yixi_time_series(hc, time, 'Julian Day', 'Top 200-m Heat Cont. (J)')
xlim([51.8271 105.7820]) % time(49:459)

xlabel('')
set(gca,'XTickLabel',[])

load qnet
hc_interped=interp1(time,hc,plot_day);

plot_day_daily=min(plot_day):max(plot_day);
for i=1:length(plot_day_daily)
    hc_daily_mean(i)=nanmean(hc_interped(floor(plot_day)==plot_day_daily(i)));
end
% hc_daily_movmean=movmean(hc_interped,8,'omitnan');
% plot(plot_day,hc_daily_movmean,'linewidth',1.7)


plot(plot_day_daily,hc_daily_mean,'linewidth',1.7)
plot(plot_day,hc_interped)

%% wind speed plot
figure
UU=met.U;
Yixi_time_series(UU, met.day, 'Year Day', 'Wind Speed (m/s)')
xlim([51.8271 105.7820]) % time(49:459)

% xlabel('')
% set(gca,'XTickLabel',[])

UU_daily_movmean=movmean(UU,8,'omitnan');
plot(met.day,UU_daily_movmean,'linewidth',1.7)


%% time series plot
figure
subplot(3,1,1)
plot_output_setup(OutputFile,'PWP_met_ERA5',ct_f,'\Theta above Freezing(ºC)')
caxis([-0.05 2.5]) 
ylim([0 600])
xlabel('')

subplot(3,1,2)
plot_output_setup(OutputFile,'PWP_met_ERA5',sa,'S_A (g kg^{-1})')
caxis([33.5 34.55])
ylim([0 600])
xlabel('')

% figure
subplot(3,1,3)
rho0=gsw_rho(sa,ct,0)-1000;
plot_output_setup(OutputFile,'PWP_met_ERA5',rho0,'\rho_0 (kg m^{-3})')
caxis([26.65 27.73]) 
ylim([0 600])

% subplot(2,1,2)
% [N2, p_mid]=gsw_Nsquared(pwp_output.s,pwp_output.t,pres,nanmean(profile.lat));
% plot_output_setup(OutputFile,'PWP_met_ERA5',[zeros(1,30);real(log10(N2))],'ln_{N^2} (rad^2 s^-^2)')
% % caxis([-18 -10])
% caxis([-14 -8])
% ylim([0 600])
