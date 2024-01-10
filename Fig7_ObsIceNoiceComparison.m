bye
%% %%%%%%%%%%%%%%%%% define region
cd ../WORKSPACE
load 2D_data_all_interp

juld=datenum(yr,mon,d)-datenum(2014,1,1)+1;

clearvars flag_depth_day_grad east west north south sal temp yr

cri4=nansum(ct')==0;
pres(isnan(ct))=NaN;


% cri=juld(:,1)<47.0038 | juld(:,1)>113.5598;
cri=juld(:,1)<50 | juld(:,1)>200;
cri1=lat(:,1)>-73 | long(:,1)>-106;
cri2=long(:,1)<-114;
% cri1=long(:,1)<-110 | long(:,1)>-108;
% cri2=lat(:,1)>-73.5 | lat(:,1)<-74;
cri=cri | cri1;
cri=cri | cri2;


ct(cri,:)=[];
sa(cri,:)=[];
pres(cri,:)=[];
lat(cri,:)=[];
long(cri,:)=[];
juld(cri,:)=[];
mon(cri,:)=[];
% sealNo(cri,:)=[];

clearvars -except ct sa pre mon long lat pres juld lat_AMS long_AMS sealNo sealName
ct_f=ct-gsw_CT_freezing(sa,pres);

load bathymetry_Amundsen

north=-73;south=-74.5;west=-115;east=-107;
cd ../EM959

Yixi_map_scatter(bathymetry,long(:,1),lat(:,1),mon(:,1),'Year Days',north,south,west,east)

% m_scatter(long(:,1),lat(:,1),23,mon(:,1),'filled')

Yixi_profile(ct',pres','ct','')
scatter(ct(:),pres(:),11,mon(:),'filled')
