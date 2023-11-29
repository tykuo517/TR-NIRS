
clc; clear;close all;

mua_ub=[0.46 0.45 0.042 0.5]; % 1/cm
mua_lb=[0.46 0.05 0.042 0.05]; % 1/cm
mus_ub=[57 64 23 350]; % 1/cm
mus_lb=[57 64 23 50]; % 1/cm

%% Random choose
num_random=7;
num_layer=4;

random_index_arr=rand(num_random,num_layer);

mua_param_arr=zeros(size(random_index_arr));
for L=1:num_layer
    mua_param_arr(:,L)=random_index_arr(:,L)*(mua_ub(L)-mua_lb(L))+mua_lb(L);
end

mus_param_arr=zeros(size(random_index_arr));
for L=1:num_layer
    mus_param_arr(:,L)=random_index_arr(:,L)*(mus_ub(L)-mus_lb(L))+mus_lb(L);
end

param_arr = [mua_param_arr mus_param_arr];

save('target_OP_to_sim.txt','param_arr','-ascii','-tabs');

%% Choose in same step

% varying mua scalp
for i=1:10
    mu(i,:) = [0.1+0.02*i 0.31 0.042 0.32 57 64 23 245];
end

% varying mua skull
for i=1:10
    mu(end+1,:) = [0.46 0.1+0.02*i 0.042 0.32 57 64 23 245];
end

% varying mua CSF
for i=1:10
    mu(end+1,:) = [0.46 0.31 0.02+0.006*i 0.32 57 64 23 245];
end

% varying mua gray matter
for i=1:10
    mu(end+1,:) = [0.46 0.31 0.042 0.1+0.02*i 57 64 23 245];
end

% varying mus gray matter
for i=1:20
    mu(end+1,:) = [0.46 0.1 0.042 0.32 57 64 23 50+i*15];
end
save('OP_to_sim.txt','mu','-ascii','-tabs');

