
% close all;

sbj_name={'KB','BY','CS','TY','CT','ZJ','WH','WW','YF','BT'};

figure('Units','pixels','Position',[0 0 1920 1080]);
ti=tiledlayout(2,5);
for i=1:length(sbj_name)
    load(fullfile('models',['headModel' sbj_name{i} '_EEG.mat']));
    Fp2=round(EEG.Fp2);
    Fp2(3)=Fp2(3)+round(1.5*10/0.93); % above Fp2 0.75 cm
    
    vol(Fp2(1),Fp2(2),Fp2(3))=7;
    
    Fp2(2)=Fp2(2)-round(1.22*10/0.93);
%     Fp2(1)=Fp2(1)-round(0.89*10/0.93);
    vol(Fp2(1),Fp2(2),Fp2(3))=7;
    
    nexttile;
    h=heatmap(squeeze(vol(:,:,Fp2(3)))');
    h.YDisplayData = flipud(h.YDisplayData);
    grid off;
    ax = gca;
    ax.XDisplayLabels = nan(size(ax.XDisplayData));
    ax.YDisplayLabels = nan(size(ax.YDisplayData));
    
    
    title(sbj_name{i});
end
title(ti,'above Fp2 1.5 cm');
% title(ti,'on Fp2');

if ~exist('results','file')
    mkdir('results');
end

print(fullfile('results','place_probe_Fp2_above_1.5cm.png'),'-dpng','-r200');
% print(fullfile('results','place_probe_Fp2.png'),'-dpng','-r200');

    

