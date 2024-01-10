% this is the most updated version, 15/Jun/2023, but Yixi expands the
% figures vertically a bit by dragging them to the same height as the
% mac pro screen
bye; 

load PWP_profile_3day
load PWP_met_ERA5_hourly_albedoed.mat

n=0;
cp=4000; %J/(kg*ºC)
day_start=[51.8297323193263];
Lf=3.35e5;
rho_FreshWater=1000;
rho_a=1.275;
cpa=1005;
Cd=0.001;
Patm=101325;
ep=0.62197;
Lv=2.501e6;
% cd dt60_dz2

OutputFile_str={'PWP_output_NoFlux','PWP_output_025fluxes','PWP_output_05fluxes','PWP_output_075fluxes','PWP_output_ControlFluxes','PWP_output_125Fluxes','PWP_output_15Fluxes'};
% OutputFile_str={'PWP_output_NoFlux','PWP_output_05fluxes','PWP_output_ControlFluxes','PWP_output_15Fluxes'};

for outputFile_num=1:length(OutputFile_str);
    n=n+1;
    OutputFile=OutputFile_str{outputFile_num};
    load(OutputFile)
    ml(n,:)=pwp_output.ml;
    pres=gsw_p_from_z(-pwp_output.z,nanmean(profile.lat));
    sa=gsw_SA_from_SP(pwp_output.s,pres,mean(profile.long),mean(profile.lat));
    SSS(n)=sa(1,end);
    hi(n,:)=pwp_output.hi;
    ct=gsw_CT_from_t(sa,pwp_output.t,pres);
    ct_f=ct-gsw_CT_freezing(sa,pres);
    rho=gsw_rho(sa,ct,pres);
    hc_ly=ct_f.*rho.*cp;
    hcn(n,:)=sum(hc_ly(1:251,:));
    sc_ly=sa.*rho.*1e-3;
    scn(n,:)=sum(sc_ly(1:251,:));
end

% mPWP_phys_EM959_advection('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_adv_fulldepth_2e_6',1,1,0,-2e-6,-2e-6)
% mPWP_phys_EM959_advection('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_adv_fulldepth_4e_6',1,1,0,-4e-6,-4e-6)
% mPWP_phys_EM959_advection('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_adv_fulldepth_6e_6',1,1,0,-6e-6,-6e-6)
% mPWP_phys_EM959_advection('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_adv_fulldepth_1e_5',1,1,0,-1e-5,-1e-5)
% mPWP_phys_EM959_advection('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_adv_fulldepth_2e_5',1,1,0,-2e-5,-2e-5)
% mPWP_phys_EM959_advection('PWP_met_ERA5', 'PWP_profile_3day_EOS', 'PWP_output_1_1_0_adv_fulldepth_4e_5',1,1,0,-4e-5,-4e-5)

n=n+1;
OutputFile='PWP_output_ControlFluxes_adv_fulldepth_2e_5';
% OutputFile='PWP_output_ControlFluxes_north';
% OutputFile='PWP_output_1_1_0_onlyMomentum';

% pwp_output.ml=pwp_output.ml(1:432);
% pwp_output.hi=pwp_output.hi(1:432);


% title(OutputFile)
load(OutputFile)
ml(n,:)=pwp_output.ml;
pres=gsw_p_from_z(-pwp_output.z,nanmean(profile.lat));
sa=gsw_SA_from_SP(pwp_output.s,pres,mean(profile.long),mean(profile.lat));
SSS(n)=sa(1,end);
hi(n,:)=pwp_output.hi;
ct=gsw_CT_from_t(sa,pwp_output.t,pres);
ct_f=ct-gsw_CT_freezing(sa,pres);
rho=gsw_rho(sa,ct,pres);
hc_ly=ct_f.*rho.*cp;
hcn(n,:)=sum(hc_ly(1:251,:));

sc_ly=sa.*rho.*1e-3;
scn(n,:)=sum(sc_ly(1:251,:));

% cd ..

% figure
% plot([0:25:100],ml(:,end),'Marker','o','LineWidth',3)
% hold on
% ylim([[0 245]])
% ylabel('Mixed-layer Depth (m)','fontsize',19); % right y-axis
% set(gca,'FontSize',17);
% xticks(0:25:100)
% xticklabels({'0%sSCI';'25%sSCI';'50%sSCI';'75%sSCI';'100%sSCI'});
% grid on
% grid minor
% line([0 100],[140 140],'LineStyle',':','LineWidth',2,'color',[239 59 167]./255)

% figure('Position',[0 0 671 721])

%% 250-m heat content
% figure('Position',[0 0 1200 721])
figure('Position',[0 0 1200 1121])

h=subplot(4,1,1);
% set(h,'position',[0.1100 0.675 0.7750 0.29])
set(h,'position',[0.1100 0.760 0.7750 0.23])

hold on
pwp_output.day=pwp_output.time+day_start;

set(gca,'FontSize',19);
set(gcf,'PaperPositionMode','auto')

ylabel('Upper 250 m Heat Content (J)','fontsize',23);
% cb=colorbar;
cmap=cbrewer2('div','Spectral',15,'PCHIP');
cmap(7:8,:)=[];
cmap=cmap(1:2:end,:);
cmap_8=cbrewer2('div','Spectral',15,'PCHIP');
cmap(4,:)=cmap_8(7,:);

% cmap(4:5,:)=[];
load Yixi_c_pink
load Yixi_c_fern

% Yixi_c_fern=d;
cmap(end+1,:)=Yixi_c_fern;

% colormap(flipud(cmap))
xlabel('Year Day','fontsize',23)

for i=1:size(hcn,1)-1
    plot(pwp_output.day,hcn(i,:),'color',cmap(i,:),'LineWidth',3)
end
    plot(pwp_output.day,hcn(end,:),'color',cmap(end,:),'LineWidth',6,'LineStyle',':')

% grid minor; box on;
xlim([min(pwp_output.day) max([pwp_output.day])])
% ylim([-0.1 1.4])
grid on; set(gca,'layer','top');
grid minor; set(gca,'layer','top')
box on; set(gca,'layer','top')
xlabel(''); xticks([50:5:110]); xticklabels('');

%%%%%%%%%%%%%%%%%% plot the obs. hc
load obs_hc % from obs_HeatBudget.m
regressionCoef_ab_hc=robustfit(obs_hc_time(obs_hc_time>day_start & obs_hc_time<max([pwp_output.day])),obs_hc(obs_hc_time>day_start & obs_hc_time<max([pwp_output.day])));
% plot(obs_hc_time(obs_hc_time>day_start),obs_hc_time(obs_hc_time>day_start)*regressionCoef_ab_hc(2)+regressionCoef_ab_hc(1))
plot(obs_hc_time,obs_hc,'color','k','LineWidth',1.3)
plot(obs_hc_days+0.5,obs_hc_daily,'color','k','LineWidth',3)
%%%%%%%%%%%%%%%%%% plot the obs. hc

ylim([0.8 4] *1e8)

% %%%%%%%%%%%%%%%%%% plot the accumulated heat flux
% OutputFile='PWP_output_1_1_0';
% load(OutputFile)
% plot(pwp_output.time+day_start,cumsum(pwp_output.qnet)*3*60*60+mean(hcn(:,1)),'color',[249 94 166]/255,'LineWidth',3)
% %%%%%%%%%%%%%%%%%% plot the accumulated heat flux
% 
% %%%%%%%%%%%%%%%%%% plot the accumulated heat for sea ice formation
% iceforming=pwp_output.hi;
% ml_depth=pwp_output.ml;
% mr=-iceforming./ml_depth;
% heatIce=mr*1000*Lf.*ml_depth;
% plot(pwp_output.time+day_start,hcn(1,:)+heatIce,'LineWidth',3)
% %%%%%%%%%%%%%%%%%% plot the accumulated heat for sea ice formation


%% top-250 salt content
h=subplot(4,1,2);
% set(h,'position',[0.1100 0.375 0.7750 0.29])
set(h,'position',[0.1100 0.52 0.7750 0.23])

hold on
pwp_output.day=pwp_output.time+met.day(1);

set(gca,'FontSize',19);
set(gcf,'PaperPositionMode','auto')

ylabel('Upper 250 m Salt Content (kg)','fontsize',23);

xlabel('Year Day','fontsize',23)

for i=1:size(scn,1)-1
    plot(pwp_output.day,scn(i,:),'color',cmap(i,:),'LineWidth',3)
end
    plot(pwp_output.day,scn(end,:),'color',cmap(end,:),'LineWidth',6,'LineStyle',':')

% grid minor; box on;
xlim([min(pwp_output.day) max([pwp_output.day])])
% ylim([-0.1 1.4])
grid on; set(gca,'layer','top');
grid minor; set(gca,'layer','top')
box on; set(gca,'layer','top')
xlabel(''); xticks([50:5:110]); xticklabels('');

%%%%%%%%%%%%%%%%%% plot the obs. sc
load obs_sc % from obs_HeatBudget
plot(obs_sc_time,obs_sc,'color','k','LineWidth',1.3)
plot(obs_sc_days+0.5,obs_sc_daily,'color','k','LineWidth',3)
%%%%%%%%%%%%%%%%%% plot the obs. sc
ylim([8778 8873])
%% sea ice thickness
h=subplot(4,1,3);
% set(h,'position',[0.1100 0.075 0.7750 0.29])
set(h,'position',[0.1100 0.28 0.7750 0.23])

set(gca,'FontSize',19);
set(gcf,'PaperPositionMode','auto')

ylabel('Ice Thickness (m)','fontsize',23);
% cb=colorbar;

xlabel('Year Day','fontsize',23)
hold on

for i=1:size(hi,1)-1
    plot(pwp_output.day,hi(i,:),'color',cmap(i,:),'LineWidth',3)
end
    plot(pwp_output.day,hi(end,:),'color',cmap(end,:),'LineWidth',6,'LineStyle',':')

% grid minor; box on;
xlim([min(pwp_output.day) max([pwp_output.day])])
ylim([-0.7 2.1])
grid on; set(gca,'layer','top');
grid minor; set(gca,'layer','top')
box on; set(gca,'layer','top')
xlabel(''); xticks([50:5:110]); xticklabels('');

%%%%%%%%%%%%%%%%%%% calcualte the extra heat that has been used to form sea ice
OutputFile='PWP_output_ControlFluxes';

load(OutputFile)
% pwp_output.day=pwp_output.time+day_start;
% 
% AccumHeat_ERA5=cumsum(pwp_output.qnet)*3*60*60+mean(hcn(:,1));

load EM959_binned_calib_3m
SST=nanmean((ct_b(:,1:2))');
SST=interp1(time_b,SST,met.day);
lw=[(0.97.*(5.67e-8)).*((SST+273.15).^4)] - 0.97.*met.lw; 
%Calculate heat fluxes
qsens=(rho_a.*cpa.*Cd).*met.U.*(SST-met.tair);
%Using WMO (2008) definition for ew
ew=6.112.*exp((17.62*(SST))./(243.12+SST)); %in hPa
ew=ew*100;
%Using AOMIP definition of specific humidity. Can't find description of where 0.378 comes from
shum_sat=(ep.*ew)./(Patm -(0.378.*ew));
qlat=((rho_a*Lv*Cd).*met.U.*(shum_sat - met.shum));
qloss=lw+qsens+qlat;
qnet=-qloss+met.sw;
%%%%%%%%%%% quick and dirty fix of soem NaNs
for i=442:446
    qnet(i)=qnet(441)+((qnet(447)-qnet(441))/6)*(i-441);
end
for i=1278:1281
    qnet(i)=qnet(1277)+((qnet(1282)-qnet(1277))/5)*(i-1277);
end

qnet(1278:1281)=qnet(148);
qnet(427)=(qnet(426)+qnet(428))/2;
%%%%%%%%%%% quick and dirty fix of soem NaNs

AccumHeat_ERA5=cumsum(qnet)*60*60+mean(hcn(:,1));
% AccumHeat_ERA5=AccumHeat_ERA5';

ResidualHeat_obs=AccumHeat_ERA5-interp1(obs_hc_time,obs_hc,met.day);

%%%%%%%%%%% quick and dirty test of half qnet's effect on the sea ice formation
% PWP_output_1_1_05=load('PWP_output_1_1_05');
% qnet=PWP_output_1_1_05.pwp_output.qnet';
% AccumHeat_ERA5=cumsum(qnet)*3*60*60+mean(hcn(:,1));
% pwp_output.day=pwp_output.time+met.day(1);
% % ResidualHeat_obs=AccumHeat_ERA5-interp1(pwp_output.day,hcn(3,:),met.day);
% ResidualHeat_obs=AccumHeat_ERA5-hcn(3,:)';
% ResidualHeat_obs=interp1(pwp_output.day,ResidualHeat_obs,met.day);
%%%%%%%%%%% quick and dirty test of half qnet's effect on the sea ice formation

IceForming_t=-ResidualHeat_obs/(Lf*rho_FreshWater);

IceForming_days=floor(min(met.day)):ceil(max(met.day));
for i=1:length(IceForming_days)
    IceForming_daily_t(i)=nanmean(IceForming_t(floor(met.day)==IceForming_days(i)));
end
plot(met.day,IceForming_t-IceForming_daily_t(2),'color','k','LineWidth',1.3)
plot(IceForming_days(2:end)+0.5,IceForming_daily_t(2:end)-IceForming_daily_t(2),'color','k','LineWidth',3)
%%%%%%%%%%%%%%%%%%% calcualte the extra heat that was used to form sea ice

%%%%%%%%%%%%%%%%%%% calcualte the extra salt that has been rejected to form sea ice
load obs_Ice_s % from obs_IceFormation_Salt_250.m

% plot(obs_Ice_s_time,obs_Ice_s-obs_Ice_s_daily(2),'color',[.5 .5 .5])
% plot(obs_Ice_s_days(2:end),obs_Ice_s_daily(2:end)-obs_Ice_s_daily(2),'color',[.5 .5 .5],'LineWidth',3)

plot(obs_Ice_s_time,obs_Ice_s-obs_Ice_s_daily(2),'color',[.61 .61 .61],'LineWidth',1.3)
plot(obs_Ice_s_days(2:end)+0.5,obs_Ice_s_daily(2:end)-obs_Ice_s_daily(2),'color',[.61 .61 .61],'LineWidth',3)

% plot(obs_Ice_s_time,obs_Ice_s,'color',[.61 .61 .61],'LineWidth',1.3)
% plot(obs_Ice_s_days(2:end)+0.5,obs_Ice_s_daily(2:end),'color',[.61 .61 .61],'LineWidth',3)
%%%%%%%%%%%%%%%%%%% calcualte the extra salt that has been rejected to form sea ice

%%%%%%%%%%%%%%%%% get the regression fit
regressionCoef_ab_saltIce=robustfit(obs_Ice_s_time,obs_Ice_s-obs_Ice_s_daily(2));
regressionCoef_ab_heatIce=robustfit(met.day,IceForming_t-IceForming_daily_t(2));
plot(obs_Ice_s_time,regressionCoef_ab_saltIce(2)*obs_Ice_s_time+regressionCoef_ab_saltIce(1),'lineWidth',1)
plot(met.day,regressionCoef_ab_heatIce(2)*met.day+regressionCoef_ab_heatIce(1),'lineWidth',1)

% regressionCoef_aa_saltIce=robustfit(obs_Ice_s_time-obs_Ice_s_time(1),obs_Ice_s-obs_Ice_s_daily(2),'','','off');
% regressionCoef_aa_heatIce=robustfit(met.day-met.day(1),IceForming_t-IceForming_daily_t(2),'','','off');
% plot(obs_Ice_s_time,regressionCoef_aa_saltIce*(obs_Ice_s_time-obs_Ice_s_time(1)),'lineWidth',1)
% plot(met.day,regressionCoef_aa_saltIce*(met.day-met.day(1)),'lineWidth',1)

%%%%%%%%%%%%%%%%%% get the regression fit

%% ml evolve
h=subplot(4,1,4);
% set(h,'position',[0.1100 0.075 0.7750 0.29])
set(h,'position',[0.1100 0.04 0.7750 0.23])

pwp_output.day=pwp_output.time+day_start;

set(gca,'FontSize',19);
set(gcf,'PaperPositionMode','auto')

ylabel('ML depth','fontsize',23);
% cb=colorbar;
% colormap(flipud(cmap))
xlabel('Year Day','fontsize',23)
hold on

for i=1:size(hcn,1)-1
    plot(pwp_output.day,ml(i,:),'color',cmap(i,:),'LineWidth',3)
end
    plot(pwp_output.day,ml(end,:),'color',cmap(end,:),'LineWidth',6,'LineStyle',':')

load obs_MLD
plot(obs_MLD_time,obs_MLD,'color','k','LineWidth',1.3)
plot(obs_MLD_days+0.5,obs_MLD_daily,'color','k','LineWidth',3)
xticks([50:5:110]);

% grid minor; box on;
xlim([min(pwp_output.day) max([pwp_output.day])])
ylim([-1 370])
grid on; set(gca,'layer','top');
grid minor; set(gca,'layer','top')
box on; set(gca,'layer','top')

