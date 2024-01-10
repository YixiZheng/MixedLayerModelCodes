bye
%% read data
% filename='2003-2010_ave.ll.icepro.data';
filename03='03.ll.monthly.icepro.antarctic.data';
filename04='04.ll.monthly.icepro.antarctic.data';
% % filename='05.ll.monthly.icepro.antarctic.data';
% % filename='06.ll.monthly.icepro.antarctic.data';
% % filename='07.ll.monthly.icepro.antarctic.data';
% filename='10.ll.monthly.icepro.antarctic.data';
%% read all the data into workspace and reshape them
fids = fopen(filename03,'r','n');

x = fread(fids,'float32');
x=reshape(x,[3600 501 1 1]); x(x==0 | x>=9.9E+32)=NaN;

SIP03=squeeze(nanmean(x(:,:,:,1:1),4));

fids = fopen(filename04,'r','n');

x = fread(fids,'float32');
x=reshape(x,[3600 501 1 1]); x(x==0 | x>=9.9E+32)=NaN;

SIP04=squeeze(nanmean(x(:,:,:,1:1),4));

SIP=(SIP03+SIP04)./2;

long=[0.:0.1:359.9]-360;
lat=-80.0:0.05:-55;


figure('Position',[771 671 871 571])
hold on

load bathymetry_Amundsen
north=-73.5;south=-74.5;west=-114;east=-108;

lat_lim=lat>=south & lat<=north;

long_lim=long>=west & long<=east;

long=long(long_lim);
lat=lat(lat_lim);
SIP=SIP(long_lim,lat_lim);

m_proj('lambert','lat',[south north],'long',[west east]);
m_pcolor(long,lat,SIP')

m_grid('fontsize',17,'tickdir','out','backgroundcolor','none');

Rtopomask=double(ncread('RTopo-2.0.1_30sec_Antarctica_aux.nc','amask'));
Rtopolong=ncread('RTopo-2.0.1_30sec_Antarctica_aux.nc','lon');
Rtopolat=ncread('RTopo-2.0.1_30sec_Antarctica_aux.nc','lat');
Rtopomask=Rtopomask(Rtopolong>=west & Rtopolong<=east,Rtopolat>=south & Rtopolat<=north);
lat_bathy=bathymetry.lat(bathymetry.lat>=south & bathymetry.lat<=north);
long_bathy=bathymetry.long(bathymetry.long>=west & bathymetry.long<=east);
m_contour(long_bathy,lat_bathy,Rtopomask',[0.99 1.99],'linecolor','k','linewidth',.3);


load EM959_binned_calib
spelim=long_b<-111 & lat_b>-74.2;
long_b=long_b(spelim);
lat_b=lat_b(spelim);
time_b=time_b(spelim);
dis_lim=[1 29 145];
long_b(dis_lim)=[];
lat_b(dis_lim)=[];
time_b(dis_lim)=[];

%% centre 
longcentre=min(long_b)+diff([min(long_b) max(long_b)])/2;
latcentre=min(lat_b)+diff([min(lat_b) max(lat_b)])/2;
%% median distance
for i=1:length(long)
    for j=1:length(lat)
        dis(i,j)=gsw_distance([long(i) longcentre],[lat(j) latcentre]);
        if dis(i,j)<=1.9*1e4
            m_scatter(long(i),lat(j),13,'Marker','*','MarkerEdgeColor','k')
        end
    end
end

colorbar

% nanmedian(dis)
load Yixi_c_pink
m_range_ring(longcentre,latcentre,19,'color',Yixi_c_pink,'LineWidth',3)

SIP_circle=SIP(dis<=1.9*1e4);

nanmean(SIP_circle)/30.5



