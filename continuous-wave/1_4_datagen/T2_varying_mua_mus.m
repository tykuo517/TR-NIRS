%{
Plot the DTOF by varying mua, mus in +-50%

Ting-Yi Kuo
Last update: 2023/6/9
%}

% load('/home/md703/Documents/Ty/1_4_MCX_lookup_table_dataGen/KB_test1_2023-05-18-16-12-12/all_param_arr.mat')
load('/home/md703/Documents/Ty/1_4_MCX_lookup_table_dataGen/KB_test1_2023-06-26-23-05-44/all_param_arr.mat')

SDS_dist_arr=[0.8 1.5 2.12 3]; % cm


%% Same mua,mus, different SDS
figure('Position',[0 0 1800 600]);
ti=tiledlayout(1,4);

for s=1:4
    nexttile;
    for i=1:10
        semilogy(1:1:10,all_param_arr(i,9+10*(s-1):18+10*(s-1)));
        hold on
    end
    title(['SDS ' num2str(s)]);
    hold off
end
title(ti,'Varying \mu_{a,scalp} for step 0.02 cm^{-1}');
xlabel('Time gate');
ylabel('Reflectance');


figure('Position',[0 0 1800 600]);
ti=tiledlayout(1,4);

for s=1:4
    nexttile;
    for i=11:20
        semilogy(1:1:10,all_param_arr(i,9+10*(s-1):18+10*(s-1)));
        hold on
    end
    title(['SDS ' num2str(s)]);
    hold off
end
title(ti,'Varying \mu_{a,skull} for step 0.02 cm^{-1}');
xlabel('Time gate');
ylabel('Reflectance');


figure('Position',[0 0 1800 600]);
ti=tiledlayout(1,4);

for s=1:4
    nexttile;
    for i=21:30
        semilogy(1:1:10,all_param_arr(i,9+10*(s-1):18+10*(s-1)));
        hold on
    end
    title(['SDS ' num2str(s)]);
    hold off
end
title(ti,'Varying \mu_{a,CSF} for step 0.006 cm^{-1}');
xlabel('Time gate');
ylabel('Reflectance');


figure('Position',[0 0 1800 600]);
ti=tiledlayout(1,4);

for s=1:4
    nexttile;
    for i=31:40
        semilogy(1:1:10,all_param_arr(i,9+10*(s-1):18+10*(s-1)));
        hold on
    end
    title(['SDS ' num2str(s)]);
    hold off
end
title(ti,'Varying \mu_{a,GM} for step 0.02 cm^{-1}');
xlabel('Time gate');
ylabel('Reflectance');

close all;

figure('Position',[0 0 1800 600]);
ti=tiledlayout(1,4);

for s=1:4
    nexttile;
    for i=41:50
        semilogy(1:1:10,all_param_arr(i,9+10*(s-1):18+10*(s-1)));
        hold on
    end
    title(['SDS ' num2str(SDS_dist_arr(s)) ' cm']);
    hold off
end
title(ti,'Varying \mu_{s,GM} for step 15 cm^{-1}');
xlabel('Time gate');
ylabel('Reflectance');

print('sensitivity_mus_GM.png','-dpng','-r200');

figure('Position',[0 0 1800 600]);
ti=tiledlayout(1,3);

%% Same SDS, different mua,mus
for s=4
    nexttile;
    for i=11:20
        semilogy(1:1:10,all_param_arr(i,9+10*(s-1):18+10*(s-1)));
        hold on
    end
    title('Varying \mu_{a,skull} for step 0.02 cm^{-1}');
    xlabel('Time gate');
    ylabel('Reflectance');
    
    nexttile;
    for i=31:40
        semilogy(1:1:10,all_param_arr(i,9+10*(s-1):18+10*(s-1)));
        hold on
    end
    title('Varying \mu_{a,GM} for step 0.02 cm^{-1}');
    xlabel('Time gate');
    ylabel('Reflectance');
    
    nexttile;
    for i=41:50
        semilogy(1:1:10,all_param_arr(i,9+10*(s-1):18+10*(s-1)));
        hold on
    end
    title('Varying \mu_{s,GM} for step 15 cm^{-1}');
    xlabel('Time gate');
    ylabel('Reflectance');

end

print('sensitivity_mus_GM.png','-dpng','-r200');




%%
% load('/home/md703/Documents/Ty/1_4_MCX_lookup_table_dataGen/KB_2023-04-28-23-58-57/all_param_arr.mat')
% 
% same_mus_index=sum(round(all_param_arr(:,5:8),3)==[57 64 23 245],2)==4;
% same_mus=all_param_arr(same_mus_index,9:18);
% 
% figure('Position',[0 0 1920 1080]);
% ti=tiledlayout(1,3);
% nexttile;
% semilogy(1:1:10,same_mus);
% title('Varying \mu_{a}');
% xlabel('Time gate');
% ylabel('Reflectance');
% 
% 
% same_mua_index=(sum(round(all_param_arr(:,1:4),3)==[0.460,0.192,0.0420,0.212],2)==4) + (mod(all_param_arr(:,8)-50,15)==0) == 2;
% % same_mua_index=sum(round(all_param_arr(:,1:4),3)==[0.460,0.192,0.0420,0.212],2)==4;
% same_mua=all_param_arr(same_mua_index,9:18);
% 
% nexttile;
% semilogy(1:1:10,same_mua);
% title('Varying \mu_{s}');
% xlabel('Time gate');
% ylabel('Reflectance');


