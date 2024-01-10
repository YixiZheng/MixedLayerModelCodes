% This is the most updated version, 12 Aug, 2023
% Yixi specially wrote this for EM959

% The data are a bit noisy, so we read the first three dayas data
% then average them to a single profile
% Yixi also movmean the profile, and make sure that the density conversion
% is gone, yet the stratification remains

% This one uses practical salinity and potential temperature

bye

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
time=repmat(time,size(pres,1),1);


last_grid=nanmean(pt(first_3,:))';
last_grid=max(find(isnan(last_grid)==0));

t_original=nanmean(pt(first_3,1:last_grid))';
s_original=nanmean(sp(first_3,1:last_grid))';
z_original=nanmean(dep(first_3,1:last_grid))';

% get a "near-surface" layer
t_original(1:2)=t_original(3);
s_original(1:2)=s_original(3);

d_original=sw_dens0(s_original,t_original);

movmean_thickness=1;
% profile.t=movmean(t_original,movmean_thickness);
% profile.s=movmean(s_original,movmean_thickness);
profile.t=t_original;
profile.s=s_original;

profile.z=-z_original;

pres_profile=nanmean(pres(first_3,1:last_grid));
profile.d=sw_dens0(profile.s,profile.t);

profile.long=repmat(nanmean(long(first_3)),last_grid,1);
profile.lat=repmat(nanmean(lat(first_3)),last_grid,1);

% % get a "over-surface" layer
% profile.t=[profile.t(1);profile.t'];
% profile.s=[profile.s(1);profile.s'];
% profile.z=[profile.z(1);-profile.z']; 

% profile.long=[profile.long(1);-profile.long];
% profile.lat=[profile.lat(1);-profile.lat];


% quick test
find(diff(profile.d)<0)
load Yixi_c_pink
load Yixi_c_blue
load Yixi_c_orange
load Yixi_c_orange_light
load Yixi_c_blue_light

ax1 = axes;
% Yixi_profile(NaN(size(profile.d)),profile.z,'\rho (kg m^-^3)','')
% ylim([0 620]);xlim([1026 1031])

ax2 = axes;
Yixi_profile(profile.d,profile.z,'\rho (kg m^-^3)','')
plot(profile.d,profile.z,'LineWidth',5,'color',Yixi_c_pink)
plot(d_original,-nanmean(dep(first_3,1:last_grid)),'k')
ylim([0 620]);xlim([1026 1031])
ax = gca;
ax.XColor = Yixi_c_pink;

ax3 = axes;
Yixi_profile(profile.t,profile.z,'Potential Temperature (ºC)','')
scatter(reshape(pt(first_3,1:last_grid),1,[]),-reshape(dep(first_3,1:last_grid),1,[]),5,'filled','MarkerFaceColor',Yixi_c_blue_light)
plot(profile.t,profile.z,'LineWidth',5,'color',Yixi_c_blue)
plot(t_original,-nanmean(dep(first_3,1:last_grid)),'k')
ylim([0 620]);xlim([-2 -0.5])
xticks(-2:0.3:-0.5)
ax = gca;
ax.XColor = Yixi_c_blue;

ax4 = axes;
Yixi_profile(profile.s,profile.z,'Salinity','')
scatter(reshape(sp(first_3,1:last_grid),1,[]),-reshape(dep(first_3,1:last_grid),1,[]),5,'filled','MarkerFaceColor',Yixi_c_orange_light)
plot(profile.s,profile.z,'LineWidth',5,'color',Yixi_c_orange)
plot(s_original,-nanmean(dep(first_3,1:last_grid)),'k')
ylim([0 620]);xlim([33 34.5])
xticks(33:0.3:34.5)
ax = gca;
ax.XColor = Yixi_c_orange;
%% link the two figures together
% linkaxes([ax1,ax2,ax3])

ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];
ax4.Visible = 'off';
ax4.XTick = [];
ax4.YTick = [];
box on
grid on
grid minor
yticks(1/620*[0:100:600])
set(gca,'YDir','reverse')


% figure
% Yixi_CT_SA(profile.t,profile.s,profile.d,profile.long,profile.lat,13,.7)
% xlim([33.1 34.75])
% ylim([-1.8 0])

save PWP_profile_3day_EOS profile

grid on
grid minor
ylim([-620 0])
box on
