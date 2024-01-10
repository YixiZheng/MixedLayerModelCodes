goodbye

cd ../../PWP/Amundsen_ERA5

%% long and lat
ERA5.long=ncread('Amundsen_hourly_all_2014.nc','longitude');
% ERA5.long=ERA5.long;
ERA5.lat=ncread('Amundsen_hourly_all_2014.nc','latitude');

%% take the 3-hourly data first (ERA5_2014_fc_sfc_Amundsen)
ERA5.time=ncread('Amundsen_hourly_all_2014.nc','time'); % start from 1.125 day
ERA5.day=(double(ERA5.time)-(datenum(2014,1,1)-datenum(1900,1,1))*24)./24+1;

ERA5.u=ncread('Amundsen_hourly_all_2014.nc','u10');
ERA5.v=ncread('Amundsen_hourly_all_2014.nc','v10');

ERA5.tair=ncread('Amundsen_hourly_all_2014.nc','t2m');
ERA5.tair=ERA5.tair-273.15;

ERA5.q=ncread('Amundsen_hourly_q_2014.nc','q');

% ERA_5.al=ncread('AMS_ERA_2014_all_south.nc','al',[1,1,1],[inf,inf,length(ERA_5.time)]);

%%
ERA5.U=sqrt(ERA5.u.^2+ERA5.v.^2);

Cd=0.001;
rho_air= 1.275;

%rho_air=1.23;
%Cd=1.3*10^(-3);
ERA5.tx=rho_air.*Cd.*abs(ERA5.u).*ERA5.u;
ERA5.ty=rho_air.*Cd.*abs(ERA5.v).*ERA5.v;

clearvars -except ERA5

%% read averaged data, ERA5_2014_fc_sfc_Amundsen_mean
namei='Amundsen_hourly_all_2014.nc';

ERA5.mean_downsw=ncread(namei,'msdwswrf');
% mean_time=ncread(namei,'time',[3],[length(time)],[3]);
ERA5.mean_downlw=ncread(namei,'msdwlwrf');

% ERA5.mean_upsh=ncread(namei,'sshf')./(60*60);
% % mean_time=ncread(namei,'time',[3],[length(time)],[3]);
% ERA5.mean_uplh=ncread(namei,'slhf')./(60*60);

ERA5.mean_precip=ncread(namei,'mtpr')./1000; % original unit is kg m**-2 s**-1, so needs to over the freshwater density, to be m/s


%% select location: DIS
% long range: [-112.125 -111.125]
% lat range: [-74.125 -73.875]
% year-day range: [51 106]

% longlim=56;
longlim=53:56;
latlim=17;
daylim=1201:2521;
ERA5.day=ERA5.day(daylim);
ERA5.time=ERA5.time(daylim);
ERA5.long=ERA5.long(longlim);
ERA5.lat=ERA5.lat(latlim);
ERA5.u=ERA5.u(longlim,latlim,daylim);
ERA5.v=ERA5.v(longlim,latlim,daylim);
ERA5.tair=ERA5.tair(longlim,latlim,daylim);
ERA5.q=ERA5.q(longlim,latlim,daylim);
ERA5.U=ERA5.U(longlim,latlim,daylim);
ERA5.tx=ERA5.tx(longlim,latlim,daylim);
ERA5.ty=ERA5.ty(longlim,latlim,daylim);
ERA5.mean_downsw=ERA5.mean_downsw(longlim,latlim,daylim);
ERA5.mean_downlw=ERA5.mean_downlw(longlim,latlim,daylim);
% ERA5.mean_upsh=ERA5.mean_upsh(longlim,latlim,daylim);
% ERA5.mean_uplh=ERA5.mean_uplh(longlim,latlim,daylim);
ERA5.mean_precip=ERA5.mean_precip(longlim,latlim,daylim);

mERA5.u=squeeze(mean(mean(ERA5.u)));
mERA5.v=squeeze(mean(mean(ERA5.v)));
mERA5.tair=squeeze(mean(mean(ERA5.tair)));
mERA5.q=squeeze(mean(mean(ERA5.q)));
mERA5.U=squeeze(mean(mean(ERA5.U)));
mERA5.tx=squeeze(mean(mean(ERA5.tx)));
mERA5.ty=squeeze(mean(mean(ERA5.ty)));
mERA5.mean_downsw=squeeze(mean(mean(ERA5.mean_downsw)));
mERA5.mean_downlw=squeeze(mean(mean(ERA5.mean_downlw)));
% mERA5.mean_upsh=squeeze(mean(mean(ERA5.mean_upsh)));
% mERA5.mean_uplh=squeeze(mean(mean(ERA5.mean_uplh)));
mERA5.mean_precip=squeeze(mean(mean(ERA5.mean_precip)));

cd ../../iSTAR/EM959
save ERA5_Dotson_hourly ERA5 mERA5



