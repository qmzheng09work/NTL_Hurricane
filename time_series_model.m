%% prepation for analysis

%实现功能
% 1. input hurricane info data:
%自动提取开始日期，hurricane_id,相关的state_id, declare disaster county_id
% 2. 模块：getpts 检视NTL pixel time series
% 3. smooth time series模块。分别smooth estimated/actual time series
% 4. 评价indicator模块：
% 1. 恢复时间
% 2. 速率
% 3. 如何判定恢复？
% 4. 如何统计恢复结果

clc;
close all;
clear;

warning('off')

path_NTL_folder='E:\OneDrive\Research\NTL_recovery_inequality\data\raster\VNP46A2_GF\';
NTL_img_list=dir([path_NTL_folder,'*.tif']);

path_urban_folder='E:\OneDrive\Research\NTL_recovery_inequality\data\raster\urban_county\';
urban_image_list=dir([path_urban_folder,'*.tif']);

path_hurricane_list='E:\OneDrive\Research\NTL_recovery_inequality\data\Hurricanes\H_Batch\';

hurricane_list=dir([path_hurricane_list,'*.txt']);
path_output='E:\OneDrive\Research\NTL_recovery_inequality\data\raster\Test_results_v3\';



for h=1:length(hurricane_list)

    Hurricane=hurricane_list(h).name(1:end-4);

    path_FEMA_record=[path_hurricane_list,Hurricane,'.txt'];
    FEMA_record=readtable(path_FEMA_record,'ReadVariableNames',1);
    
    state_list=unique(FEMA_record.State);

    for s=1:length(state_list)
        target_state=state_list{s};

        record_county=[];
        para=[];

        path_output_folder=[path_output,Hurricane,'_',target_state,'\'];
        if exist(path_output_folder)==0
            mkdir(path_output_folder)
        end

        idx=zeros(1,size(FEMA_record,1))';

        for i=1:length(idx)
            if strcmp(cell2mat(FEMA_record{i,6}),target_state)
                idx(i)=1;
            end
        end

        FEMA_target_state=FEMA_record(find(idx>0),:);

        state_id=FEMA_target_state{1,7};
        state_id=floor(state_id/100);
        dec_county_id_list=FEMA_target_state{:,7};

        t1 = datetime(2012,1,20);

        time_str=num2str(unique(FEMA_record.date));
        t_event=datetime(str2double(time_str(1:4)),str2double(time_str(5:6)),str2double(time_str(7:8)));
        t_event_idx=split(between(t1,t_event,'Days'),'days');

        record_county.Hurricane=Hurricane;
        record_county.target_state=target_state;
        record_county.dec_county_id_list=dec_county_id_list;
        record_county.FEMA_type=FEMA_target_state{:,8};


        thre=30;
        time_buffer=30;
        prob_thre=0.95;

        para.thre=thre;
        para.time_buffer=time_buffer;
        para.prob_thre=prob_thre;

        record_county.para=para;

        for c=1:length(dec_county_id_list)

            county_id=dec_county_id_list(c);

            NTL_img=double(imread([path_NTL_folder,'VNP46A2_GF_',num2str(county_id),'.tif']));
            [urban_img,R]=geotiffread([path_urban_folder,'urban_',num2str(county_id),'.tif']);
            urban_img=double(urban_img==13);


            [m,n]=size(urban_img);

            NTL_img_lit_area=double(median(NTL_img(:,:,1:t_event_idx),3,'omitnan')>thre);
            NTL_img_thre=NTL_img;

            parfor i=1:size(NTL_img,3)
                NTL_img_thre(:,:,i)=NTL_img(:,:,i).*NTL_img_lit_area;
            end

            NTL_ts_sum=zeros(1,size(NTL_img,3));
            for i=1:size(NTL_img,3)
                NTL_ts_sum(i)=sum(NTL_img_thre(:,:,i),'all');
            end

            try
                BEAST_output=beast_irreg(medfilt1(NTL_ts_sum(:),7),'start', 1, 'season','harmonic','freq',365,'tseg.min',300,'mcmc.seed',3,'print.progress',0,'print.options',0);
                [flag_county,match_timing_county,recovery_time_county,trend_cp_prob_county]=fCheck_Change_Occurrence_v1(BEAST_output,t_event_idx,time_buffer,prob_thre);

                record_county.flag_county=flag_county;
                record_county.recovery_time_county=recovery_time_county;
                record_county.trend_cp_prob_county=trend_cp_prob_county;
            catch
                error_idx_county=1;
                record_county.error_idx_county=error_idx_county;
            end

            change_img=zeros(m,n);
            trend_cp_prob_img=zeros(m,n);
            recovery_duration=zeros(m,n);
            error_idx=zeros(m,n);

            for i=1:m
                parfor j=1:n
                    if NTL_img_lit_area(i,j)==0
                        continue
                    end

                    NTL_ts_selected=NTL_img_thre(i,j,:);

                    try
                        BEAST_output=beast_irreg(medfilt1(NTL_ts_selected(:),7), 'start', 1, 'season','harmonic','freq',365,'tseg.min',300,'mcmc.seed',3,'print.progress',0,'print.options',0);
                        [flag,match_timing,recovery_time,trend_cp_prob]=fCheck_Change_Occurrence_v1(BEAST_output,t_event_idx,time_buffer,prob_thre);

                        if flag==0
                            change_img(i,j)=1; 
                        elseif flag==-1
                            change_img(i,j)=3; 
                        else
                            change_img(i,j)=2; 
                            recovery_duration(i,j)=recovery_time;
                            trend_cp_prob_img(i,j)=trend_cp_prob;
                        end
                    catch
                        error_idx(i,j)=1;
                    end
                end
            end

            record_county.county_id=county_id;
            record_county.t_event=t_event;
            record_county.t_event_idx=t_event_idx;
            record_county.change_img=change_img;
            record_county_trend_cp_prob_img=trend_cp_prob_img;
            record_county.per_change_pixel=length(find(change_img==2))/length(find(change_img>0));

            record_county.R=R;
            record_county.recovery_duration=recovery_duration;

            num_change_pix=length(find(change_img==2));
            num_long_change_pixel=length(find(recovery_duration<0))+length(find(recovery_duration>recovery_time_county));
            record_county.per_longer_recovery_pix=num_long_change_pixel/num_change_pix;
            record_county.error_idx=error_idx;



            % save and output
            save([path_output_folder,Hurricane,'_',num2str(county_id),'_record.mat'],"record_county");
            geotiffwrite([path_output_folder,Hurricane,'_',num2str(county_id),'_change_img.tif'],change_img,R);
            geotiffwrite([path_output_folder,Hurricane,'_',num2str(county_id),'_recovery_duration.tif'],recovery_duration,R);

        end
    end
end