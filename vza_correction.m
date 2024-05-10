clear;
close all;
clc;

path_root='E:\OneDrive\Research\NTL_recovery_inequality\data\raster\';

path_BRDF=[path_root,'VNP46A2_BRDF\'];
path_VZA=[path_root,'VNP46A2_VZA\'];
path_GF=[path_root,'VNP46A2_GF\'];

path_output=[path_root,'VNP46A2_GF_corrected\'];

state_id_all=[12,13,22,37,45,48,72];
thre_NTL=15; % 10 nW/sr/cm2
thresh_quantile=0.15;

record_all=[];

for s=1:length(state_id_all)
    state_id=state_id_all(s);

    img_list=dir([path_VZA,'*',num2str(state_id),'*.tif']);
    for c=1:length(img_list)


        VZA_img_name=img_list(c).name;
        GF_img_name=strrep(VZA_img_name,'VZA','GF');
        BRDF_img_name=strrep(VZA_img_name,'VZA','BRDF');
        GF_img_corrected_name=strrep(VZA_img_name,'VZA','GF_cor');
        county_id=VZA_img_name(end-8:end-4);

        record_all(c).state_id=c;
        record_all(c).county_id=county_id;

        [VZA_img,R]=geotiffread([path_VZA,VZA_img_name]);
        VZA_img=double(VZA_img);
        GF_img=double(imread([path_GF,GF_img_name]));
        BRDF_img=double(imread([path_BRDF,BRDF_img_name]));

        [row,col,T]=size(GF_img);

        t_train=1:1000;
        t_validation=(length(t_train)+1):T;
        GF_img_corrected=GF_img;
        avg_rad=mean(GF_img,3,'omitnan');

        record_D_metric=zeros(row*col,2);
        record_stat=zeros(row*col,2);
        record_delta_NTL=zeros(row*col,length(0:90));

        for i =1:(row*col)

            [x,y]=ind2sub([row,col],i);
            if avg_rad(x,y)<=thre_NTL
                continue
            end

            BRDF_ts=BRDF_img(x,y,t_train);
            GF_ts=GF_img(x,y,t_train);
            VZA_ts=VZA_img(x,y,t_train);
            [b,stats,delta_NTL]=fVZA_calibration(GF_ts,BRDF_ts,VZA_ts,thresh_quantile);

            record_delta_NTL(i,:)=delta_NTL;
            record_stat(i,1)=stats(1);
            record_stat(i,2)=stats(3);


            [D_metric,D_metric_cal,GF_ts_corrected]=fVZA_correction(b,x,y,BRDF_img,GF_img,VZA_img,0:70);
            record_D_metric(i,1)=D_metric;
            record_D_metric(i,2)=D_metric_cal;

            GF_img_corrected(x,y,:)=GF_ts_corrected;
        end

        record_all(c).stat=record_stat;
        record_all(c).D_metric=record_D_metric;
        record_all(c).delta_NTL=record_delta_NTL;

        geotiffwrite([path_output,GF_img_corrected_name],GF_img_corrected,R);

    end
end


