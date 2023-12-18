%{
Using the pre-simulated ANN to forward spectrum, and compare the result to the target
only consider the error of choosed SDS

Ting-Yi Kuo
Last update: 2023/08/21
%}

function output=fun_forward_calError_chooseSDS(param_arr)
    global testing_index output_folder target_spec orig_target_spec target_dtof final fitRoute_txt fitRoute_csv times;
    global lambda cw_flag tr_flag fitting_wl_tr num_SDS_cw num_SDS_tr SDS_choosed_cw SDS_choosed_tr SDS_sim_correspond_cw SDS_sim_correspond_tr num_gate gate_choosed gate_sim_correspond use_min_disk mask weight;  

    %% Forward the spec
    [OP_arr,~]=fun_param_to_mu(param_arr,0);

    if cw_flag
        ann_spec=fun_ANN_forward(OP_arr,0);
        ann_spec=ann_spec(:,SDS_sim_correspond_cw); % only use the SDS that are in target spectrum
    end
    
    if tr_flag
        OP_arr_tr=interp1(lambda,OP_arr,fitting_wl_tr,'pchip');
        ann_dtof_=fun_ANN_forward(OP_arr_tr,1);
        
        for i=1:num_SDS_tr
            ann_dtof(:,i)=ann_dtof_((i-1)*num_gate+1:(i-1)*num_gate+num_gate)';
        end
        ann_dtof=ann_dtof(:,SDS_sim_correspond_tr); % only use the SDS that are in target spectrum
    end
    
    %% Calculate error
    if cw_flag
        SDS_error=sqrt(mean((ann_spec./target_spec-1).^2,1));
        error(1)=sqrt(mean(SDS_error(SDS_choosed_cw).^2));
        error(1)=error(1)*weight(1)/sum(weight);
    else
        SDS_error=NaN*ones(1,num_SDS_cw);
        error(1)=NaN;
    end
    
    if tr_flag
        select_gate=sum(mask,1);
        temp_SDS_error=sqrt((ann_dtof./target_dtof-1).^2);
        temp_SDS_error=temp_SDS_error.*mask;
        temp_SDS_error=sqrt(sum(temp_SDS_error.^2,1)./select_gate);

        SDS_error(end+1:end+num_SDS_tr)=temp_SDS_error;
        error(2)=sqrt(mean(temp_SDS_error(SDS_choosed_tr).^2));
        error(2)=error(2)*weight(2)/sum(weight);
    else
        SDS_error(end+1:end+num_SDS_tr)=NaN*ones(1,num_SDS_tr);
        error(2)=NaN;
    end
    
    % total error
    if sum(isnan(error))==0
        choosed_error=sum(error);
    else
        choosed_error=sum(error(~isnan(error)));
    end
    SDS_error(end+1:end+3)=[error choosed_error];
    
    SDS_error_disp=[SDS_error(1:num_SDS_cw+num_SDS_tr)];
    
    % disp error
    fprintf('Error for each SDS:');
    fprintf('\t%.2f%%',SDS_error_disp*100);
    fprintf(',\tmean = %.2f%%\n',SDS_error(num_SDS_cw+num_SDS_tr+3)*100);
    
    
    %% Save the fitting route
    if testing_index==1
        fitRoute_txt = fopen(fullfile(output_folder,'fitRoute.txt'), 'a');
        fitRoute_csv = fopen(fullfile(output_folder,'fitRoute.csv'), 'a');
        while fitRoute_txt==-1
            fitRoute_txt = fopen(fullfile(output_folder,'fitRoute.txt'), 'a');
        end
        while fitRoute_csv==-1
            fitRoute_csv = fopen(fullfile(output_folder,'fitRoute.csv'), 'a');
        end
        fprintf(fitRoute_csv,'A1,K1,A2,K2,A4,K4,hc_1,StO2_1,hc_2,StO2_2,hc_4,StO2_4,mel,');
        for s=1:num_SDS_cw
            fprintf(fitRoute_csv,'error %d,',s);
        end
        for s=1:num_SDS_tr
            fprintf(fitRoute_csv,'error %d,',s);
        end
        fprintf(fitRoute_csv,'cw error,tr error,mean error,\n');
    end
    
    if use_min_disk==0 || testing_index==1 || final==1
        if testing_index>1
            fprintf(fitRoute_txt,'\n');
            fprintf(fitRoute_csv,'\n');
        end
        
        param_error_set = [param_arr SDS_error];
        fprintf(fitRoute_txt,'%f ', param_error_set);
        fprintf(fitRoute_csv,'%f,', param_error_set);
    end
    
    if final==0
        testing_index=testing_index+1;
        output=choosed_error;
    else
        fclose(fitRoute_csv);
        fclose(fitRoute_txt);
        % plot and save the spec
        if use_min_disk==0
            fig=figure('Units','pixels','position',[0 0 1920 1080]);
            set(fig, 'visible', 'off');
            ti=tiledlayout('flow','TileSpacing','compact','Padding','none'); 
            for s=1:num_SDS
                nexttile();
                hold on;
                plot(orig_target_spec(:,1),orig_target_spec(:,s+1));
                plot(lambda,ann_spec(:,s));
                title(['SDS ' num2str(s) ', error = ' num2str(SDS_error(s)*100,'%.2f%%')]);
            end
            title(ti,['try ' num2str(testing_index)]);
            saveas(gcf,fullfile(output_folder,'fitted_spec.png'));
            close all;
        end
        
        if cw_flag
            save(fullfile(output_folder,'fitted_spec.txt'),'ann_spec','-ascii','-tabs');
        end
        if tr_flag
            save(fullfile(output_folder,'fitted_dtof.txt'),'ann_dtof','-ascii','-tabs');
        end
        
        save(fullfile(output_folder,'fitted_mu.txt'),'OP_arr','-ascii','-tabs');
        fprintf('Total fitting iteration: %d\n',times);
    end
end

%% functions
