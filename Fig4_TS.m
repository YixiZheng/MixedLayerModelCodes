goodbye
load EM959_calib
% scatter(sa,ct,3,'filled','MarkerFaceColor',[.7 .7 .7])
% hold on
% clear

% spelim=long<-111 & (lat>-74.15 & lat<-73.87);
spelim=long<-111 & lat>-74.2;
long=long(spelim);
lat=lat(spelim);
time=time(spelim);
ct=ct(spelim);
sa=sa(spelim);
dep=dep(spelim);

timelim=time>70.2 & time<70.35;
long(timelim)=[];
lat(timelim)=[];
time(timelim)=[];
ct(timelim)=[];
sa(timelim)=[];
dep(timelim)=[];
%% TS plot for QC
% figure('Position',[0 0 471 421])
figure('Position',[0 0 571 521])

Yixi_CT_SA(ct,sa,-dep,17,.7)
xlim([33.22 34.85])
ylim([-2 .7])
% xlim([33.2 34.8]); ylim([-2 .4])
caxis([0 600])
% scatter3(sa,time,ct,13,dep,'filled');
% colormap(cbrewer('div','Spectral',16,'PCHIP'))
% colorbar
% xlim([33.2 34.8]); zlim([-2 .4])

figure('Position',[0 0 571 521])
Yixi_CT_SA(ct,sa,time,17,.7)
% xlim([32.7 35])
xlim([33.22 34.85])
ylim([-2 .7])
caxis([min(time) 133])% the date of the last winter ones
% 

longlim=[min(long) max(long)];
latlim=[min(lat) max(lat)];

clearvars -except longlim latlim
cd ../WORKSPACE
load 2D_data_AMS
cd ../EM959

%% test plot
% load bathymetry_AMS
% north=-73;south=-75;west=-115;east=-100;

% figure
% Yixi_map_scatter(bathymetry,long(mon(:,1)>4,:),lat(mon(:,1)>4,:),mon(mon(:,1)>4,:),'month',north,south,west,east)

%% more winter dots
timelim=(mon(:,1)>=4 & sealNo(:,1)~=5)& (sum(isnan(sa'))<18)';
% timelim=mon(:,1)>=4 & (sum(isnan(sa'))<18)';
spalim=(long(:,1)>=longlim(1) & long(:,1)<=longlim(2)) & (lat(:,1)>=latlim(1) & lat(:,1)<=latlim(2));
sain=sa(timelim & spalim,:);
ctin=ct(timelim & spalim,:);
monin=mon(timelim & spalim,:);
sealNoin=sealNo(timelim & spalim,:);

dep=gsw_z_from_p(pres,lat);
depin=dep(timelim & spalim,:);

days=datenum(yr,mon,d)-datenum(2014,1,1)+1;
daysin=days(timelim & spalim,:);

sc=scatter(sain(monin==5),ctin(monin==5),13,daysin(monin==5));
plot((sain(monin(:,1)==9,:))',(ctin(monin(:,1)==9,:))','Color',[.3 .3 .3]);
scatter(sain(monin==9),ctin(monin==9),17,'filled','MarkerFaceColor',[.3 .3 .3]);

colormap(cbrewer2('div','Spectral',51,'PCHIP'))
colorbar('northoutside','color','k')
colorbar('eastoutside','color','k')

set(gca,'color',[.8 .8 .8])
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'color','w')

figure(1)
sc=scatter(sain(monin==5),ctin(monin==5),13,-depin(monin==5));
plot((sain(monin(:,1)==9,:))',(ctin(monin(:,1)==9,:))','Color',[.3 .3 .3]);
scatter(sain(monin==9),ctin(monin==9),17,-depin(monin==9));
colormap(cbrewer2('div','Spectral',51,'PCHIP'))
colorbar('northoutside')

set(gca,'color',[.8 .8 .8])
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'color','w')

