%{
Try to generate noise with slope, and plot them

Benjamin Kao
Last update: 2021/01/21
%}

clc;clear;close all;

%% param
num_error_to_generate=100; % number of adding noise to the same, the first one will have no error
num_error_to_plot=10;
SDS_error_arr=[3 4.2 5.1 5.2 5.4 12.1]; % in %, the error (CV) of each SDS
small_noise_ratio=0.05;
num_SDS=6;
SDS_dist_arr=[0.8 1.5 2.12 3 3.35 4.5 4.74]; % cm
wl=(700:900)';

fontSize=12;

%% main

error_arr=[];
for i=1:num_error_to_generate
    for s=1:num_SDS
        % decide the 2 end of the slope line
        error_end=normrnd(0,SDS_error_arr(s),1,2);
        error_end(2)=error_end(1)+(error_end(2)-error_end(1))*rand(1);
        SDS_error_spec=interp1([min(wl) max(wl)],error_end,wl);

        % add small noise
        SDS_error_spec=SDS_error_spec+normrnd(0,SDS_error_arr(s)*small_noise_ratio,size(wl));
        error_arr(:,s,i)=SDS_error_spec;
    end
end

fig=figure('Units','pixels','position',[0 0 1200 800]);
ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
for s=1:num_SDS
    nexttile();
    hold on;
    for i=1:num_error_to_plot
        plot(wl,error_arr(:,s,i));
    end
    xlabel('wavelength(nm)');
    ylabel('error(%)');
    grid on;
    title(['SDS = ' num2str(SDS_dist_arr(s)) ' cm, mean CV= ' num2str(mean(std(error_arr(:,s,:),[],3)),'%.2f') '%']);
    
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
end

print('added_Slope_Error.png','-dpng','-r200');