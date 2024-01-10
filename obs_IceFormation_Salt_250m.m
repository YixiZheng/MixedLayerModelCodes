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
% hc=heat_content(sa_hc,ct_hc,pres,1,interp_depth);

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
interp_depth=250;

ct=ct(spelim,1:interp_depth+1);
sa=sa(spelim,1:interp_depth+1);
time=time(spelim);
long=long(spelim);
lat=lat(spelim);
pres=pres(spelim,1:interp_depth+1);

sa(sum((isnan(ct(:,1:interp_depth+1)))')>1,:)=[];
pres(sum((isnan(ct(:,1:interp_depth+1)))')>1,:)=[];
time(sum((isnan(ct(:,1:interp_depth+1)))')>1)=[];
ct(sum((isnan(ct(:,1:interp_depth+1)))')>1,:)=[];
%%%%%%%%%%%%%%%%%%% if they don't reach 250 m, omit
ct_f=gsw_CT_freezing(sa,pres);
% rho=gsw_rho(sa,ct,pres);
rho=gsw_rho(sa,ct,pres);

sarho=(mean((sa.*rho)'))';

%%%%%%%%%%%%%%%%%%% calculate salt content
% sc=(sum((sa./1e3.*rho)'))';
% 
% plot_day_daily=floor(min(time)):ceil(max(time));
% for i=1:length(plot_day_daily)
%     sc_daily_mean(i)=nanmean(sc(floor(time)==plot_day_daily(i)));
% end
% plot(plot_day_daily,sc_daily_mean,'linewidth',1.7)
% obs_sc=sc;obs_sc_time=time;
% obs_sc_daily=sc_daily_mean'; obs_sc_days=plot_day_daily';
% obs_sc_sa=sa;
% clearvars -except obs_sc obs_sc_daily obs_sc_days obs_sc_time obs_sc_sa
% save obs_sc
%%%%%%%%%%%%%%%%%%% calculate salt content


for i=2:size(ct,1)
    Vi(i)=(interp_depth*(sarho(i)-sarho(i-1)))./(sarho(i)-920*7);
end

obs_ice_test=(interp_depth*(sarho-sarho(1)))./(sarho-920*7);% just to make sure that the following code is correct

obs_Ice_s=cumsum(Vi);obs_Ice_s_time=time;

obs_Ice_s_time(obs_Ice_s==NaN)=[];
obs_Ice_s(obs_Ice_s==NaN)=[];

% subplot(3,1,1)
Yixi_time_series(obs_Ice_s, obs_Ice_s_time, 'Year Day', 'Sea Ice Formation')

obs_Ice_s_days=floor(min(time)):ceil(max(time));
for i=1:length(obs_Ice_s_days)
    obs_Ice_s_daily(i)=nanmean(obs_Ice_s(floor(time)==obs_Ice_s_days(i)));
end
plot(obs_Ice_s_days,obs_Ice_s_daily,'linewidth',1.7)
legend('daily data','7-day Running Average')

clearvars -except obs_Ice_s_days obs_Ice_s_daily obs_Ice_s obs_Ice_s_time
save obs_Ice_s

