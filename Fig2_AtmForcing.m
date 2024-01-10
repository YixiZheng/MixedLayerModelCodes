bye
% format long
load Yixi_c_PrettyFive.mat
load Yixi_c_ColourPack.mat
load PWP_profile_3day
load PWP_met_ERA5_hourly_albedoed.mat
load timelim

% 
% %%%%%%%%%%%%%%%%%%% calcualte the extra heat that has been used to form sea ice
% iceforming=pwp_output.hi;
% ml_depth=pwp_output.ml;
% mr=-iceforming./ml_depth;
Lf=3.35e5;
% Heatice=mr*1000*Lf.*ml_depth;
% qice=diff(Heatice)/(3*60*60);
Lv=2.501e6;
% 
% n=0;
rho_FreshWater=1000;
rho_Ice=920;
rho_a=1.275;
cpa=1005;
Cd=0.001;
Patm=101325;
ep=0.62197;
% %%%%%%%%%%%%%%%%%%% calcualte the extra heat that has been used to form sea ice

% %%%%%%%%%%%%%%%%%%% calcualte the heat content and heat content change
% pres=gsw_p_from_z(-pwp_output.z,nanmean(profile.lat));
% sa=gsw_SA_from_SP(pwp_output.s,pres,mean(profile.long),mean(profile.lat));
% ct=gsw_CT_from_t(sa,pwp_output.t,pres);
% ct_f=ct-gsw_CT_freezing(sa,pres);
% rho=gsw_rho(sa,ct,pres);
cp=4183.3; %J/(kg*ºC)
% hc_ly=ct_f.*rho.*cp;
% hc=sum(hc_ly)*2;
% qhc=diff(hc)/(3*60*60);
% %%%%%%%%%%%%%%%%%%% calcualte the heat content and heat content change


%%%%%%%%%%%%%%%%%%% calculate the latent and sensible heat fluxes
load EM959_binned_calib_3m
SST=nanmean((ct_b(:,1:2))');
SST=interp1(time_b,SST,met.day);
lw=[(0.97.*(5.67e-8)).*((SST+273.15).^4)] - 0.97.*met.lw; 
qsens=(rho_a.*cpa.*Cd).*met.U.*(SST-met.tair);
%Using WMO (2008) definition for ew
ew=6.112.*exp((17.62*(SST))./(243.12+SST)); %in hPa
ew=ew*100;
%Using AOMIP definition of specific humidity. Can't find description of where 0.378 comes from
shum_sat=(ep.*ew)./(Patm -(0.378.*ew));
qlat=((rho_a*Lv*Cd).*met.U.*(shum_sat - met.shum));
qloss=lw+qsens+qlat;
qnet=-qloss+met.sw;
%%%%%%%%%%%%%%%%%%% calculate the latent and sensible heat fluxes

% figure('Position',[0 0 1071 1021])
% figure('Position',[0 0 1071 821])
figure('Position',[0 0 1071 821]*1.2)

%% latent, sensible, sw, lw
% h=subplot(5,1,1);
% set(h,'position',[0.1100 0.82 0.7750 0.17])
h=subplot(3,1,1);
set(h,'position',[0.1100 2/3+0.02-0.01 0.7750 0.315])
hold on

set(gca,'FontSize',19);
set(gcf,'PaperPositionMode','auto')

ylabel('Heat Flux (W m^-^2)','fontsize',23);
xlabel('Year Day','fontsize',23)

% plot(met.day,-met.sw,'LineWidth',1.3)
% plot(met.day,-met.lw,'LineWidth',1.3,'color',Yixi_c_pink)
% ylim([-450 1])

line([min(met.day) max([met.day])],[0 0],'color','k','LineStyle',':','LineWidth',2)
plot(met.day,qnet,'LineWidth',3,'color',[.67 .67 .67])
plot(met.day,met.sw,'LineWidth',1.3,'color',Yixi_c_PrettyFive(1,:))
plot(met.day,lw,'LineWidth',1.3,'color',Yixi_c_PrettyFive(3,:))
plot(met.day,-qlat,'LineWidth',1.3,'color',Yixi_c_PrettyFive(4,:))
plot(met.day,-qsens,'LineWidth',1.3,'color',Yixi_c_PrettyFive(5,:))
ylim([-390 450])

xlim([min(met.day) max([met.day])])
grid on; grid minor
box on; 
xlabel(''); xticks([50:5:110]); xticklabels('');

% legend('Shortwave Radiation (W m^-^2)','Longwave Radiation (W m^-^2)')

% %% 2-m air temp
% h=subplot(5,1,2);
% set(h,'position',[0.1100 0.64 0.7750 0.17])
% hold on
% 
% set(gca,'FontSize',19);
% set(gcf,'PaperPositionMode','auto')
% 
% xlabel('Year Day','fontsize',23)
% 
% plot(met.day,met.tair,'LineWidth',1.3)
% ylim([-18 -.1])
% ylabel('Air Temperature (ºC)','fontsize',23,'lineWidth',1.3)
% 
% xlim([min(met.day) max([met.day])])
% grid on; grid minor
% box on; 
% xticks([50:5:110]); 
% xticklabels(''); 
% 
% % ylim([-400 300])

%% wind stress 
% h=subplot(5,1,3);
% set(h,'position',[0.1100 0.46 0.7750 0.17])
h=subplot(3,1,2);
set(h,'position',[0.1100 1/3+0.03-0.01 0.7750 0.315])
hold on

set(gca,'FontSize',19);
set(gcf,'PaperPositionMode','auto')

xlabel('Year Day','fontsize',23)
% 
% xfact=48;
% yfact=1;
% quiver([met.day;55;55],[zeros(size(met.day));-0.05;-0.05],[met.tx*xfact;.2*xfact;0],[met.ty*yfact;0;.2*yfact],0,'ShowArrowHead','off','LineWidth',1,'AutoScale','off')
% text(55,-0.1,'0.2 m s^-^2','FontSize',19)

% yyaxis left
% ylim([-0.3 0.3])
% plot(met.day,met.tx,'LineWidth',1)
% 
% yyaxis right
ylim([-0.31 0.299])
% plot(met.day,met.ty,'LineWidth',1)


plot(met.day,met.tx,'LineWidth',1.3)
plot(met.day,met.ty,'LineWidth',1.3,'color',Yixi_c_pink)
plot(met.day,met.U./100,'LineWidth',2.9,'color',Yixi_c_orange)

line([50 110],[0 0],'color','k','LineWidth',2,'LineStyle',':')
% legend('Eastward Wind Stress (m s^-^2)','Northward Wind Stress (m s^-^2)')
ylabel('Wind Stress (m s^-^2)','fontsize',23)

ylabel('Wind Stress and Speed','fontsize',23)

xlim([min(met.day) max([met.day])])
grid on; grid minor
box on; 
xticks([50:5:110]); 
xticklabels(''); 

% ylim([-400 300])

%% P-E
% h=subplot(5,1,4);
% set(h,'position',[0.1100 0.28 0.7750 0.17])
h=subplot(3,1,3);
set(h,'position',[0.1100 0/3+0.04-0.01 0.7750 0.315])

hold on
set(gca,'FontSize',19);
set(gcf,'PaperPositionMode','auto')

line([min(met.day) max([met.day])],[0 0],'color','k','LineStyle',':','LineWidth',2)
plot(met.day,met.precip*1e7,'LineWidth',1.3,'color',Yixi_c_blue)
ylabel('Freshwater Flux (10^-^7 m s^-^1)','fontsize',23)


%Calculate evaporation and emp
evap=(qlat/(1000*(Lv)))*1e7;
emp = met.precip*1e7-evap;

plot(met.day,emp,'LineWidth',3.7,'color',[.6 .6 .6])
plot(met.day,-evap,'LineWidth',1.3,'color',Yixi_c_pink)

xlim([min(met.day) max([met.day])])
grid on; grid minor
box on; 
xticks([50:5:110]); 
ylim([-0.5 4.88])
% xticklabels(''); 

% ylim([-4.88 0.1])

% %% humidity
% h=subplot(5,1,5);
% set(h,'position',[0.1100 0.1 0.7750 0.17])
% 
% hold on
% set(gca,'FontSize',19);
% set(gcf,'PaperPositionMode','auto')
% 
% plot(met.day,met.shum*1e3,'LineWidth',1.3)
% ylabel('Specific humidity (g kg^-^1)','fontsize',23)
% 
% xlim([min(met.day) max([met.day])])
% grid on; grid minor
% box on; 
% xticks([50:5:110]); 
% % ylim([-4.88 0.1])



