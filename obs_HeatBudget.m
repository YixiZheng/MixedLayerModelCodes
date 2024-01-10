% by far the most updated one, 18/Oct/2023
bye; clf

load EM959_binned_calib_m
ct=ct_b;sa=sa_b;long=long_b;lat=lat_b;pres=pres_b;time=time_b;dep=dep_b;
clearvars ct_b sa_b long_b lat_b pres_b dep_b time_b
spelim=long<-111 & lat>-74.2;
time_bgn=min(time(spelim));
time_end=max(time(spelim));

% ct_hc=ct;sa_hc=sa;
% for i=1:size(ct,1)
%     for j=2:size(ct,2)
%         if isnan(ct_hc(i,j))==1 && isnan(ct(i,j-1))==0
%             ct_hc(i,j)=mean([ct(i,j-1),ct(i,j+1)]);
%             sa_hc(i,j)=mean([sa(i,j-1),sa(i,j+1)]);
%         end
%     end
% end
% hc=heat_content(sa_hc,ct_hc,pres,1,200);

%%%%%%%%%%%%%%%%%%% if sa<10, omit
ct(sa<10)=NaN;
rho(sa<10)=NaN;
sa(sa<10)=NaN;
%%%%%%%%%%%%%%%%%%% if sa<10, omit

%%%%%%%%%%%%%%%%%%% interp blank measurement, vertically
% fill single gaps (1 NaN)
for i=1:length(time)
    for j=2:size(sa,2)-1
        if isnan(sa(i,j))==1
            if isnan(sa(i,j-1))==0 && isnan(sa(i,j+1))==0
                sa(i,j)=0.5*(sa(i,j-1)+sa(i,j+1));
            end
            if isnan(ct(i,j-1))==0 && isnan(ct(i,j+1))==0
                ct(i,j)=0.5*(ct(i,j-1)+ct(i,j+1));
            end
        end
    end
end

% fill double gaps (2 NaNs together)
for i=1:length(time)
    for j=2:size(sa,2)-2
        if isnan(sa(i,j))==1
            if isnan(sa(i,j+1))==0
                sa(i,j)=0.5*(sa(i,j-1)+sa(i,j+1));
            elseif isnan(sa(i,j+1))==1 && isnan(sa(i,j+2))==0
                sa(i,j)=0.5*(sa(i,j-1)+sa(i,j+2));
                sa(i,j+1)=0.5*(sa(i,j-1)+sa(i,j+2));
            end

            if isnan(ct(i,j+1))==0
                ct(i,j)=0.5*(ct(i,j-1)+ct(i,j+1));
            elseif isnan(ct(i,j+1))==1 && isnan(ct(i,j+2))==0
                ct(i,j)=0.5*(ct(i,j-1)+ct(i,j+2));
                ct(i,j+1)=0.5*(ct(i,j-1)+ct(i,j+2));
            end
        end
    end
end

% fill everything to surface
for i=1:length(time)
    if min(find(isnan(ct(i,:))==0))<=5
    ct(i,1:min(find(isnan(ct(i,:))==0))-1)=ct(i,min(find(isnan(ct(i,:))==0)));
    end

    if min(find(isnan(sa(i,:))==0))<=5
    sa(i,1:min(find(isnan(sa(i,:))==0))-1)=sa(i,min(find(isnan(sa(i,:))==0)));
    end
end
%%%%%%%%%%%%%%%%%%% interp blank measurement, vertically

%%%%%%%%%%%%%%%%%%% if they don't reach 250 m, omit
sa(sum((isnan(ct(:,1:251)))')>1,:)=[];
pres(sum((isnan(ct(:,1:251)))')>1,:)=[];
time(sum((isnan(ct(:,1:251)))')>1)=[];
ct(sum((isnan(ct(:,1:251)))')>1,:)=[];
%%%%%%%%%%%%%%%%%%% if they don't reach 250 m, omit

%%%%%%%%%%%%%%%%%%% calculate hc
ct_f=gsw_CT_freezing(sa,pres);
rho=gsw_rho(sa,ct,pres);
cp=4000; %J/(kg*ºC)
hc_pres=(ct(:,1:251)-ct_f(:,1:251)).*rho(:,1:251).*cp;
hc=NaN(size(pres,1),1);
hc=sum(hc_pres');
%%%%%%%%%%%%%%%%%%% calculate hc
time(isnan(hc)==1)=[];
hc(isnan(hc)==1)=[];

subplot(3,1,1)
Yixi_time_series(hc, time, 'Year Day', 'Top 200-m Heat Cont. (J)')
xlim([51.8271 105.7820]) % time(49:459)

xlabel('')
set(gca,'XTickLabel',[])

%%%%%%%%%%%%%%%% this part use the three-hourly interped hc
% load qnet
% hc_interped=interp1(time,hc,plot_day);
%
% plot_day_daily=min(plot_day):max(plot_day);
% for i=1:length(plot_day_daily)
%     hc_daily_mean(i)=nanmean(hc_interped(floor(plot_day)==plot_day_daily(i)));
% end
% hc_daily_movmean=movmean(hc_interped,8,'omitnan');
% plot(plot_day,hc_daily_movmean,'linewidth',1.7)
%%%%%%%%%%%%%%%% this part use the three-hourly interped hc

%%%%%%%%%%%%%%%% this part use the daily-bin-averaged hc
plot_day_daily=floor(min(time)):floor(max(time));
for i=1:length(plot_day_daily)
    hc_daily_mean(i)=mean(hc(floor(time)==plot_day_daily(i)));
    profile_num(:,:,i)=size(hc(floor(time)==plot_day_daily(i)));
end
profile_num=squeeze(profile_num);
median(profile_num(2,:));
%%%%%%%%%%%%%%%% this part use the daily-bin-averaged hc


%%%%%%%%%%%%%%%%%%%% saving the hc dataset for Fig.5
obs_hc=hc;obs_hc_time=time;
obs_hc_daily=hc_daily_mean'; obs_hc_days=plot_day_daily';
clearvars -except obs_hc obs_hc_daily obs_hc_days obs_hc_time
save obs_hc
%%%%%%%%%%%%%%%%%%%% saving the hc dataset for Fig.5

% plot(plot_day_daily,hc_daily_mean,'linewidth',1.7,'color',[.3 .3 .3])
% time_hc=time;
% % save hc hc time_hc

% ylim([1.2 2.7]*1e8)

%% heat flux from heat content
% % the unit of heat flux (W m^-2) has been converted to 3-hour heat input (J m^-2)
% hc_loss_daily=diff(hc_daily_mean)/(24*60*60);
% hc_loss_interped=diff(hc_interped)/(3*60*60);
% subplot(3,1,2)
% Yixi_time_series(hc_loss_interped, plot_day(1:end-1)+diff(plot_day)/2, 'Julian Day', 'Heat Flux from Seals (W m^-^2)')
% % Yixi_time_series([], [], 'Julian Day', 'Heat Flux from Seals (W m^-^2)')
% plot(plot_day_daily(1:end-1)+diff(plot_day_daily)/2,hc_loss_daily,'linewidth',1.7)
% 
% % hi_mean=movmean(hc_loss_interped,56,'omitnan');
% % plot(plot_day(1:end-1)+diff(plot_day)/2,hi_mean,'linewidth',1.7)
% % legend('daily Running Average data','7-day Running Average')
% line([51.8271 105.7820],[0 0],'linewidth',1.3,'color','k')
% xlim([51.8271 105.7820]) % time(49:459)
% ylim([-5 5]*1e3)
% xlabel('')
% set(gca,'XTickLabel',[])
% % the extra heat loss (accumulated atm. net heat loss - upper ocean net
% % heat loss) is used for forming sea ice
