% a note, this version is the most up-to-date version, 22/Dec/22, 18/Apr/2023
bye
figure('Position',[771 671 871 571])

load bathymetry_Amundsen
north=-73.5;south=-74.5;west=-114;east=-108;

bedrock=bathymetry.bedrock(bathymetry.long>=west & bathymetry.long<=east,bathymetry.lat>=south & bathymetry.lat<=north);
ice=bathymetry.ice(bathymetry.long>=west & bathymetry.long<=east,bathymetry.lat>=south & bathymetry.lat<=north);
lat_bathy=bathymetry.lat(bathymetry.lat>=south & bathymetry.lat<=north);
long_bathy=bathymetry.long(bathymetry.long>=west & bathymetry.long<=east);

for i = 1:size(bedrock,1)
    for j = 1:size(bedrock,2)
        if bedrock(i,j)>0
            bedrock(i,j)=NaN;
            ice(i,j)=NaN;
        end
    end
end

%% plot the bathymetry map
ax1 = axes;

m_proj('lambert','lat',[south north],'long',[west east]);

m_pcolor(long_bathy,lat_bathy,bedrock');
m_grid('fontsize',17,'tickdir','out','backgroundcolor','none');
caxis(ax1,[-1500 0])
colormap(ax1,gray)

%% land mask
ax2 = axes;
m_proj('lambert','lat',[south north],'long',[west east]);
Rtopomask=double(ncread('RTopo-2.0.1_30sec_Antarctica_aux.nc','amask'));
Rtopolong=ncread('RTopo-2.0.1_30sec_Antarctica_aux.nc','lon');
Rtopolat=ncread('RTopo-2.0.1_30sec_Antarctica_aux.nc','lat');
Rtopomask=Rtopomask(Rtopolong>=west & Rtopolong<=east,Rtopolat>=south & Rtopolat<=north);
m_contourf(long_bathy,lat_bathy,Rtopomask',[0.99 1.99],'linecolor','k','linewidth',.3);
colormap(ax2,[.7 .7 .7;1 1 1]);
m_grid('fontsize',17,'tickdir','out','backgroundcolor','none');

%% ice map
[fig_modis, R] = geotiffread('WorldView20140314_AMS.tif');
info = geotiffinfo('WorldView20140314_AMS.tif');

height = info.Height; % Integer indicating the height of the image in pixels
width = info.Width; % Integer indicating the width of the image in pixels
[rows,cols] = meshgrid(1:height,1:width);

% Getting the latitude and longitude of the points
[x1,y1] = pix2map(info.RefMatrix, rows, cols);
[y,x] = projinv(info,y1,x1);

x=info.CornerCoords.Lon(1)+x(1,1)-x;

load ice_map
load ice_lucky_interp
ax3 = axes;
m_proj('lambert','lat',[south north],'long',[west east])
fig_modis_ice=(double(fig_modis)).^0.3;
% fig_modis_ice(fig_modis_ice<300)=NaN;
m_pcolor(x,y,fig_modis_ice');
% m_pcolor(x,y,double(fig_modis'))
% hold on
% fig_modis_below210=fig_modis;
% fig_modis_below210(fig_modis_below210>=210)=NaN;
% m_pcolor(x,y,fig_modis_below210);
% colormap(ax3,flipud(ice_map))
colormap(ax3,gray)
caxis(ax3,[5 14])
%caxis(ax3,[3 14])

% caxis(ax3,[-2500 7000])
m_grid('fontsize',17,'tickdir','out','backgroundcolor','none');
alpha(ax3,1)

%% EM959
load EM959_binned_calib
ax4 = axes;
hold on
colormap(ax4,cbrewer2('div','Spectral',ceil(max(time_b)-min(time_b)),'PCHIP'))
caxis([min(time_b) max(time_b)])

spelim=long_b<-111 & lat_b>-74.2;
long_b=long_b(spelim);
lat_b=lat_b(spelim);
time_b=time_b(spelim);
dis_lim=[1 29 145];
long_b(dis_lim)=[];
lat_b(dis_lim)=[];
time_b(dis_lim)=[];

m_scatter(long_b,lat_b,29,time_b,'filled');

m_grid('fontsize',17,'tickdir','out','backgroundcolor','none');

m_contour(long_bathy,lat_bathy,Rtopomask',[0.99 1.99],'linecolor','k','linewidth',.3);

%% wind speed, added 18/Apr/2023
load ERA5_ASP_wind
wind_u=squeeze(mean(ERA5_ASP_wind.u(1:2:end,:,:),3,'omitnan'));
wind_v=squeeze(mean(ERA5_ASP_wind.v(1:2:end,:,:),3,'omitnan'));
[wind_long,wind_lat]=meshgrid(ERA5_ASP_wind.long(1:2:end),ERA5_ASP_wind.lat);
wind_long=wind_long';wind_lat=wind_lat';
ScaleFac=.03;
load Yixi_c_orange
m_quiver([wind_long(:);-109;-109],[wind_lat(:);-74.125;-74.125],[wind_u(:);...
    3;0]*ScaleFac,[wind_v(:);0;3]*ScaleFac,'linewidth',1.3,'color',Yixi_c_orange,'AutoScale','off');
m_quiver([-109;-109],[-74.125;-74.125],[3;0]*ScaleFac,[0;3]*ScaleFac,'linewidth',3,'color',Yixi_c_orange,'AutoScale','off');

%% link the two figures together
linkaxes([ax1,ax2,ax3,ax4])
ax3.Visible = 'off';
ax4.Visible = 'off';


ax3.XTick = [];
ax3.YTick = [];
ax4.XTick = [];
ax4.YTick = [];
set([ax1,ax2,ax3,ax4],'Position',[.225 .07 .685 .815]);
% cb1 = colorbar(ax1,'Position',[.115 .12 .02 .71],'fontsize',13,'ticks',-3000:200:0,'ticklabels',num2cell([3000:-200:-200,0]),'color','k');
cb4 = colorbar('southoutside','peer',ax4,'Position',[.23 .84 .685 .03],'fontsize',13,'color','k');

%% centre 
longcentre=min(long_b)+diff([min(long_b) max(long_b)])/2;
latcentre=min(lat_b)+diff([min(lat_b) max(lat_b)])/2;
%% median distance
for i=1:length(long_b)
    dis(i)=gsw_distance([long_b(i) longcentre],[lat_b(i) latcentre]);
end
% nanmedian(dis)
load Yixi_c_pink
m_range_ring(longcentre,latcentre,19,'color',Yixi_c_pink,'LineWidth',3)

clearvars -except time_b
cd ../WORKSPACE
load 2D_data_AMS    
cd ../EM959

%% more winter dots
longlim=[-112.1160 -111.0258];
latlim=[-74.1756 -73.859];
timelim=(mon(:,1)>=4 & sealNo(:,1)~=5);
% timelim=mon(:,1)>=4 & (sum(isnan(sa'))<18)';
spalim=(long(:,1)>=longlim(1) & long(:,1)<=longlim(2)) & (lat(:,1)>=latlim(1) & lat(:,1)<=latlim(2));
sain=sa(timelim & spalim,:);
ctin=ct(timelim & spalim,:);
longin=long(timelim & spalim,:);
latin=lat(timelim & spalim,:);
monin=mon(timelim & spalim,:);
sealNoin=sealNo(timelim & spalim,:);

days=datenum(yr,mon,d)-datenum(2014,1,1)+1;
daysin=days(timelim & spalim,:);

m_scatter(longin(monin==5),latin(monin==5),29,daysin(monin==5),'filled','MarkerEdgeColor','w');
m_scatter(longin(monin==5),latin(monin==5),29,daysin(monin==5),'MarkerEdgeColor','w');
m_scatter(longin(monin==9),latin(monin==9),29,'MarkerFaceColor',[.2 .2 .2],'MarkerEdgeColor','w','Marker','^');

caxis([min(time_b) max(daysin(monin==5))])

m_line([-112.125 -111.125],[-74.125 -74.125],'Color','w','LineWidth',1,'LineStyle','-.')
m_line([-112.125 -111.125],[-73.875 -73.875],'Color','w','LineWidth',1,'LineStyle','-.')
m_line([-112.125 -112.125],[-74.125 -73.875],'Color','w','LineWidth',1,'LineStyle','-.')
m_line([-111.125 -111.125],[-74.125 -73.875],'Color','w','LineWidth',1,'LineStyle','-.')

m_line([-111.875 -111.875],[-74.125 -73.875],'Color','w','LineWidth',1,'LineStyle','-.')
m_line([-111.625 -111.625],[-74.125 -73.875],'Color','w','LineWidth',1,'LineStyle','-.')
m_line([-111.375 -111.375],[-74.125 -73.875],'Color','w','LineWidth',1,'LineStyle','-.')

set(gcf,'color','w')

% ax = gca;
% ax.XColor = 'w';
% ax.YColor = 'w'; 
% ax.GridColor = 'k';
% ax.MinorGridColor = 'k';
% set(gcf,'color','k')
% tobj = findall(gcf,'Type','Text');
% for i = 1:length(tobj)
%     tobj(i).Color = 'w';
% end
set(gcf, 'InvertHardcopy', 'off')
m_ruler([.8 .9],.1,2,'fontsize',13);
