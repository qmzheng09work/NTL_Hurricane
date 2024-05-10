function [D_metric,D_metric_cal,GF_ts_cal]=fVZA_correction(b_est,x,y,BRDF_img,GF_img,VZA_img,time_range)

BRDF_ts_val=BRDF_img(x,y,:);
BRDF_ts_val=BRDF_ts_val(:);

GF_ts_val=GF_img(x,y,:);
GF_ts_val=GF_ts_val(:);

VZA_ts_cal=VZA_img(x,y,:);
VZA_ts_cal=VZA_ts_cal(:);
idx=double(BRDF_ts_val==GF_ts_val);
GF_ts_cal=GF_ts_val;

for v=1:length(idx)
    if idx(v) <=0
        continue
    end
    b_new=[b_est(1)/b_est(3),b_est(2)/b_est(3),0];
    GF_ts_cal(v)=BRDF_ts_val(v)/(1+polyval(b_new,VZA_ts_cal(v)));
end


VZA_ts_HighQ=VZA_ts_cal(find(idx>0));
GF_ts_HighQ=GF_ts_val(find(idx>0));
GF_ts_cal_HighQ=GF_ts_cal(find(idx>0));

X = [ VZA_ts_HighQ.^2,VZA_ts_HighQ,ones(size(VZA_ts_HighQ))];
[b,bint,r,rint,stats]=regress(GF_ts_HighQ,X);
b_new=[b(1)/b(3),b(2)/b(3),0];
D_metric=sum(abs(polyval(b_new,time_range)));

X = [ VZA_ts_HighQ.^2,VZA_ts_HighQ,ones(size(VZA_ts_HighQ))];
[b,bint,r,rint,stats]=regress(GF_ts_cal_HighQ,X);
b_cal_new=[b(1)/b(3),b(2)/b(3),0];
D_metric_cal=sum(abs(polyval(b_cal_new,time_range)));
