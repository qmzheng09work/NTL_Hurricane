function [b,stats,delta_NTL]=fVZA_calibration(GF_ts,BRDF_ts,VZA_ts,thresh)

BRDF_ts=BRDF_ts(:);
GF_ts=GF_ts(:);
VZA_ts=VZA_ts(:);

idx=double(BRDF_ts==GF_ts);

% VZA analysis

VZA_ts_HighQ=VZA_ts(find(idx>0));
BRDF_ts_HighQ=BRDF_ts(find(idx>0));

[BRDF_ts_HighQ_sort,I]=sort(BRDF_ts_HighQ);
VZA_ts_HighQ_sort=VZA_ts_HighQ(I);

num_obs=length(BRDF_ts_HighQ_sort);
BRDF_ts_HighQ_sort=BRDF_ts_HighQ_sort(floor(thresh*num_obs):floor((1-thresh)*num_obs));
VZA_ts_HighQ_sort=VZA_ts_HighQ_sort(floor(thresh*num_obs):floor((1-thresh)*num_obs));

X = [ VZA_ts_HighQ_sort.^2,VZA_ts_HighQ_sort,ones(size(VZA_ts_HighQ_sort))];
[b,bint,r,rint,stats]=regress(BRDF_ts_HighQ_sort,X);

y2=polyval(b,0:90);
y2_VZA0=y2(1);

delta_NTL=(y2-y2_VZA0)./y2_VZA0;

