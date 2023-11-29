
clear all;

global tr_net tr_param_range;
% load ANN model
load(fullfile('model_arrange',['KB_cw_model.mat'])); % cw_net, cw_param_range
load(fullfile('model_arrange',['KB_tr_model.mat'])); % tr_net, tr_param_range

change=[-50,-40,-30,-20,-10,0,10,20,30,40,50];
% change=[-20,-10,0,10,20];

baseline=[0.2,150,0.2,150,0.042,23,0.2,150];
% baseline=[0.25,150,0.2,150,0.042,23,0.25,150];

mu=[];
% test scalp mua
for i=1:length(change)
    temp=baseline(8)*(1+change(i)/100);
    mu(end+1,:)=[0.2,150,0.2,150,0.042,23,0.2,temp];
end

% mu=[0.25,154,0.36,150,0.042,23,0.36,150];
dtof=fun_ANN_forward(mu,1);

figure;
tiledlayout('flow');
for i=1:5
    nexttile;
    to_plot=dtof(:,1+10*(i-1):10*i);
    semilogy(1:1:10,to_plot);
%     ylim([min(to_plot(:)) max(to_plot(:))]);
end

relative_change=dtof./dtof(6,:)-1;
figure;
tiledlayout('flow');
for i=1:5
    nexttile;
    to_plot=relative_change(:,1+10*(i-1):6+10*(i-1));
    plot(1:1:6,to_plot);
    ylim([-0.4 0.4])
%     ylim([min(to_plot(:)) max(to_plot(:))]);
end