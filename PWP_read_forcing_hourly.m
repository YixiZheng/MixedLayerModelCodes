bye

time_bgn=51.8297323193263;
time_end=105.776610844383;
load ERA5_Dotson_hourly.mat


met.long=mean(ERA5.long);
met.lat=mean(ERA5.lat);
met.precip=squeeze(mean(mean(ERA5.mean_precip(:,:,21:end),1),2));
met.tair=squeeze(mean(mean(ERA5.tair(:,:,21:end),1),2));
met.shum=squeeze(mean(mean(ERA5.q(:,:,21:end),1),2));
met.sw=squeeze(mean(mean(ERA5.mean_downsw(:,:,21:end),1),2));
met.lw=squeeze(mean(mean(ERA5.mean_downlw(:,:,21:end),1),2));
met.day=ERA5.day(21:end);
met.time=met.day-met.day(1);
met.U=squeeze(mean(mean(ERA5.U(:,:,21:end),1),2));
met.tx=squeeze(mean(mean(ERA5.tx(:,:,21:end),1),2));
met.ty=squeeze(mean(mean(ERA5.ty(:,:,21:end),1),2));

wind.u=squeeze(mean(mean(ERA5.u(:,:,21:end),1),2));
wind.v=squeeze(mean(mean(ERA5.v(:,:,21:end),1),2));

save PWP_met_ERA5_hourly met

%% for albedo
met.sw=met.sw*0.94;% assume albedo is 0.06

save PWP_met_ERA5_hourly_albedoed met

%% for zero tau
met.tx=met.tx*0;
met.ty=met.tx*0;
save PWP_met_ERA5_hourly_0tau_albedoed met

%% for zero U
met.U=met.U.*0;

save PWP_met_ERA5_hourly_0tau_0U_albedoed met
