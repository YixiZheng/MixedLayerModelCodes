% This is the most updated version, 12 Aug, 2023
% Yixi specially wrote this for EM959

% The data are a bit noisy, so we read the first three dayas data
% then average them to a single profile
% Yixi also movmean the profile, and make sure that the density conversion
% is gone, yet the stratification remains

% This one uses practical salinity and potential temperature

goodbye

figure('position',[1340 466 1141 673])
load EM959_binned_calib_2m
time=time_b;
pres=pres_b;
long=long_b;
lat=lat_b;
dep=dep_b;
sp=gsw_SP_from_SA(sa_b,pres,long,lat);
pt=gsw_pt_from_CT(sa_b,ct_b);

clearvars time_b sa_b ct_b pres_b long_b lat_b dep_b

spelim=long<-111 & lat>-74.2;
time_bgn=min(time(spelim));
time_end=max(time(spelim));

dep=repmat(dep,size(pres,1),1);

first_3=find(time(:,1)>=time_bgn & time(:,1)<time_bgn+3);
final_10=find(time(:,1)>=time_end-3 & time(:,1)<time_end);
time=repmat(time,1,size(pres,2));


bottom_grid_first=nanmean(pt(first_3,:))'; % average potential temperature for each level
bottom_grid_first=max(find(isnan(bottom_grid_first)==0)); % the maximum level where there're still data
 
t_original=nanmean(pt(first_3,1:bottom_grid_first))';
s_original=nanmean(sp(first_3,1:bottom_grid_first))';
z_original=nanmean(dep(first_3,1:bottom_grid_first))';

bottom_grid_final=nanmean(pt(final_10,:))'; % average potential temperature for each level
bottom_grid_final=max(find(isnan(bottom_grid_final)==0)); % the maximum level where there're still data
% t_final=nanmean(pt(final_10,1:bottom_grid_final))';
% s_final=nanmean(sp(final_10,1:bottom_grid_final))';
% z_final=nanmean(dep(final_10,1:bottom_grid_final))';
t_final_mean=nanmean(pt(final_10,1:bottom_grid_final))';
s_final_mean=nanmedian(sp(final_10,1:bottom_grid_final))';
z_final_mean=nanmedian(dep(final_10,1:bottom_grid_final))';
t_final=(pt(final_10,1:bottom_grid_final))';
s_final=(sp(final_10,1:bottom_grid_final))';
z_final=(dep(final_10,1:bottom_grid_final))';

% get a "near-surface" layer
t_original(1:2)=t_original(3);
s_original(1:2)=s_original(3);

profile.t=t_original;
profile.s=s_original;

profile.z=-z_original;

pres_profile=nanmean(pres(first_3,1:bottom_grid_first));
pres_final_mean=nanmean(pres(final_10,1:bottom_grid_final));
pres_final=(pres(final_10,1:bottom_grid_final));

profile.long=repmat(nanmean(long(first_3)),bottom_grid_first,1);
profile.lat=repmat(nanmean(lat(first_3)),bottom_grid_first,1);


load Yixi_c_fern
load Yixi_c_blue
load Yixi_c_PrettyFive.mat
% cmap=cbrewer2('div','Spectral',11,'PCHIP');
% cmap(6:7,:)=[];
cmap=cbrewer2('div','Spectral',15,'PCHIP');
cmap(7:8,:)=[];
cmap=cmap(1:2:end,:);
cmap_8=cbrewer2('div','Spectral',15,'PCHIP');
cmap(4,:)=cmap_8(7,:);


PWP_output_Double=load('PWP_output_DoubleFluxes');
PWP_output_175=load('PWP_output_175Fluxes');
PWP_output_15=load('PWP_output_15Fluxes');
PWP_output_125=load('PWP_output_125Fluxes');
PWP_output_Control=load('PWP_output_ControlFluxes');
PWP_output_025=load('PWP_output_025Fluxes');
PWP_output_05=load('PWP_output_05Fluxes');
PWP_output_075=load('PWP_output_075Fluxes');
PWP_output_0=load('PWP_output_NoFlux');
PWP_output_adv=load('PWP_output_ControlFluxes_adv_fulldepth_2e_5');

profile.sa=gsw_SA_from_SP(profile.s,pres_profile',mean(profile.long),mean(profile.lat));
profile.ctf=gsw_CT_from_pt(profile.sa,profile.t)-gsw_CT_freezing(profile.sa,pres_profile');

sa_final=gsw_SA_from_SP(s_final,pres_final',mean(profile.long),mean(profile.lat));
ctf_final=gsw_CT_from_pt(sa_final,t_final)-gsw_CT_freezing(sa_final,pres_final');

sa_final_mean=gsw_SA_from_SP(s_final_mean,pres_final_mean',mean(profile.long),mean(profile.lat));
ctf_final_mean=gsw_CT_from_pt(sa_final_mean,t_final_mean)-gsw_CT_freezing(sa_final_mean,pres_final_mean');

pres=gsw_p_from_z(-PWP_output_Control.pwp_output.z,mean(profile.lat));
% PWP_output_Double.pwp_output.sa=gsw_SA_from_SP(PWP_output_Double.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
% PWP_output_Double.pwp_output.ctf=gsw_CT_from_pt(PWP_output_Double.pwp_output.sa,PWP_output_Double.pwp_output.t)-gsw_CT_freezing(PWP_output_Double.pwp_output.sa,pres);
% PWP_output_175.pwp_output.sa=gsw_SA_from_SP(PWP_output_175.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
% PWP_output_175.pwp_output.ctf=gsw_CT_from_pt(PWP_output_175.pwp_output.sa,PWP_output_175.pwp_output.t)-gsw_CT_freezing(PWP_output_175.pwp_output.sa,pres);
PWP_output_15.pwp_output.sa=gsw_SA_from_SP(PWP_output_15.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
PWP_output_15.pwp_output.ctf=gsw_CT_from_pt(PWP_output_15.pwp_output.sa,PWP_output_15.pwp_output.t)-gsw_CT_freezing(PWP_output_15.pwp_output.sa,pres);
PWP_output_125.pwp_output.sa=gsw_SA_from_SP(PWP_output_125.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
PWP_output_125.pwp_output.ctf=gsw_CT_from_pt(PWP_output_125.pwp_output.sa,PWP_output_125.pwp_output.t)-gsw_CT_freezing(PWP_output_125.pwp_output.sa,pres);
PWP_output_Control.pwp_output.sa=gsw_SA_from_SP(PWP_output_Control.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
PWP_output_Control.pwp_output.ctf=gsw_CT_from_pt(PWP_output_Control.pwp_output.sa,PWP_output_Control.pwp_output.t)-gsw_CT_freezing(PWP_output_Control.pwp_output.sa,pres);
PWP_output_025.pwp_output.sa=gsw_SA_from_SP(PWP_output_025.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
PWP_output_025.pwp_output.ctf=gsw_CT_from_pt(PWP_output_025.pwp_output.sa,PWP_output_025.pwp_output.t)-gsw_CT_freezing(PWP_output_025.pwp_output.sa,pres);
PWP_output_05.pwp_output.sa=gsw_SA_from_SP(PWP_output_05.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
PWP_output_05.pwp_output.ctf=gsw_CT_from_pt(PWP_output_05.pwp_output.sa,PWP_output_05.pwp_output.t)-gsw_CT_freezing(PWP_output_05.pwp_output.sa,pres);
PWP_output_075.pwp_output.sa=gsw_SA_from_SP(PWP_output_075.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
PWP_output_075.pwp_output.ctf=gsw_CT_from_pt(PWP_output_075.pwp_output.sa,PWP_output_075.pwp_output.t)-gsw_CT_freezing(PWP_output_075.pwp_output.sa,pres);
PWP_output_0.pwp_output.sa=gsw_SA_from_SP(PWP_output_0.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
PWP_output_0.pwp_output.ctf=gsw_CT_from_pt(PWP_output_0.pwp_output.sa,PWP_output_0.pwp_output.t)-gsw_CT_freezing(PWP_output_0.pwp_output.sa,pres);
PWP_output_adv.pwp_output.sa=gsw_SA_from_SP(PWP_output_adv.pwp_output.s,pres,mean(profile.long),mean(profile.lat));
PWP_output_adv.pwp_output.ctf=gsw_CT_from_pt(PWP_output_adv.pwp_output.sa,PWP_output_adv.pwp_output.t)-gsw_CT_freezing(PWP_output_adv.pwp_output.sa,pres);


% Yixi_profile(NaN(size(profile.d)),profile.z,'\rho (kg m^-^3)','')
% ylim([0 620]);xlim([1026 1031])

% ax2 = axes;
% Yixi_profile(profile.d,profile.z,'\rho (kg m^-^3)','')
% plot(profile.s,profile.z,'LineWidth',3,'color','k')
% plot(PWP_output_1_1_0.pwp_output.d(:,end),PWP_output_1_1_0.pwp_output.z,'LineWidth',3,'color',Yixi_c_PrettyFive(1,:))
% plot(PWP_output_1_1_025.pwp_output.d(:,end),PWP_output_1_1_025.pwp_output.z,'LineWidth',3,'color',Yixi_c_PrettyFive(2,:))
% plot(PWP_output_1_1_05.pwp_output.d(:,end),PWP_output_1_1_05.pwp_output.z,'LineWidth',3,'color',Yixi_c_PrettyFive(3,:))
% plot(PWP_output_1_1_075.pwp_output.d(:,end),PWP_output_1_1_075.pwp_output.z,'LineWidth',3,'color',Yixi_c_PrettyFive(4,:))
% plot(PWP_output_1_1_1.pwp_output.d(:,end),PWP_output_1_1_1.pwp_output.z,'LineWidth',3,'color',Yixi_c_PrettyFive(5,:))
% 
% ylim([0 620]);xlim([1026 1031])
% ax = gca;
% ax.XColor = Yixi_c_pink;
subplot(1,2,1)

Yixi_profile(profile.ctf,profile.z,'\Theta above Freezing (ºC)','')
plot(profile.ctf,profile.z,'LineWidth',3,'color','k')
plot(ctf_final,-z_final,'LineWidth',.1,'color',[.75 .75 .75])
plot(ctf_final_mean,-z_final_mean,'LineWidth',3,'color',[.6 .6 .6])
% plot(PWP_output_Double.pwp_output.ctf(:,end),PWP_output_Double.pwp_output.z,'LineWidth',3,'color',cmap(9,:),'lineStyle','-')
% plot(PWP_output_175.pwp_output.ctf(:,end),PWP_output_175.pwp_output.z,'LineWidth',3,'color',cmap(8,:),'lineStyle','-')
plot(PWP_output_15.pwp_output.ctf(:,end),PWP_output_15.pwp_output.z,'LineWidth',3,'color',cmap(7,:),'lineStyle','-')
plot(PWP_output_125.pwp_output.ctf(:,end),PWP_output_125.pwp_output.z,'LineWidth',3,'color',cmap(6,:),'lineStyle','-')
plot(PWP_output_Control.pwp_output.ctf(:,end),PWP_output_Control.pwp_output.z,'LineWidth',3,'color',cmap(5,:),'lineStyle','-')
plot(PWP_output_075.pwp_output.ctf(:,end),PWP_output_075.pwp_output.z,'LineWidth',3,'color',cmap(4,:),'lineStyle','-')
plot(PWP_output_05.pwp_output.ctf(:,end),PWP_output_05.pwp_output.z,'LineWidth',3,'color',cmap(3,:),'lineStyle','-')
plot(PWP_output_025.pwp_output.ctf(:,end),PWP_output_025.pwp_output.z,'LineWidth',3,'color',cmap(2,:),'lineStyle','-')
plot(PWP_output_0.pwp_output.ctf(:,end),PWP_output_0.pwp_output.z,'LineWidth',3,'color',cmap(1,:),'lineStyle','-')
plot(PWP_output_adv.pwp_output.ctf(:,end),PWP_output_adv.pwp_output.z,'LineWidth',5,'color',Yixi_c_fern,'lineStyle',':')

ylim([0 400]);xlim([-0.1 1.8])
% xticks(0:0.3:2)


subplot(1,2,2)
Yixi_profile(profile.sa,profile.z,'S_A (g kg^-^3)','')
plot(profile.sa,profile.z,'LineWidth',3,'color','k')
plot(sa_final,-z_final,'LineWidth',.1,'color',[.75 .75 .75])
plot(sa_final_mean,-z_final_mean,'LineWidth',3,'color',[.6 .6 .6])

% plot(PWP_output_Double.pwp_output.sa(:,end),PWP_output_Double.pwp_output.z,'LineWidth',3,'color',cmap(9,:),'lineStyle','-')
% plot(PWP_output_175.pwp_output.sa(:,end),PWP_output_175.pwp_output.z,'LineWidth',3,'color',cmap(8,:),'lineStyle','-')
plot(PWP_output_15.pwp_output.sa(:,end),PWP_output_15.pwp_output.z,'LineWidth',3,'color',cmap(7,:),'lineStyle','-')
plot(PWP_output_125.pwp_output.sa(:,end),PWP_output_125.pwp_output.z,'LineWidth',3,'color',cmap(6,:),'lineStyle','-')
plot(PWP_output_Control.pwp_output.sa(:,end),PWP_output_Control.pwp_output.z,'LineWidth',3,'color',cmap(5,:),'lineStyle','-')
plot(PWP_output_075.pwp_output.sa(:,end),PWP_output_075.pwp_output.z,'LineWidth',3,'color',cmap(4,:),'lineStyle','-')
plot(PWP_output_05.pwp_output.sa(:,end),PWP_output_05.pwp_output.z,'LineWidth',3,'color',cmap(3,:),'lineStyle','-')
plot(PWP_output_025.pwp_output.sa(:,end),PWP_output_025.pwp_output.z,'LineWidth',3,'color',cmap(2,:),'lineStyle','-')
plot(PWP_output_0.pwp_output.sa(:,end),PWP_output_0.pwp_output.z,'LineWidth',3,'color',cmap(1,:),'lineStyle','-')
plot(PWP_output_adv.pwp_output.sa(:,end),PWP_output_adv.pwp_output.z,'LineWidth',5,'color',Yixi_c_fern,'lineStyle',':')

ylim([0 400]);xlim([33.3 34.6])
% xticks(33:0.3:34.5)



