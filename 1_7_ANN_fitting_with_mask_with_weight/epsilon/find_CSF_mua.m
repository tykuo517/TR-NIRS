% water: mua of water
% hb: mua of hb
% hbo: mua of hbo

water=load('water_1_cm.txt');
hb=load('Hb_1_cm_um_mua.txt');
hbo=load('HbO_1_cm_um_mua.txt');

reference_wl=800;
refreence_mua=0.02; % cm^-1
Arachnoid_fraction=0.15;
Arachnoid_mua=0.16; % cm^-1
sto2=0.75;
alpha=0; % the fraction of hb & hbo

wl=(ceil(max([water(1,1) hb(1,1) hbo(1,1)])):floor(min([water(end,1) hb(end,1) hbo(end,1)])))';

hb=interp1(hb(:,1),hb(:,2),wl)*10^6; % turn mua(1/uM) into mua(1/M) 
hbo=interp1(hbo(:,1),hbo(:,2),wl)*10^6;

water=interp1(water(:,1),water(:,2),wl);
ua_hc = ( (sto2.*hbo/64532)+((1-sto2).*hb/64500) )*150;

index=find(wl==reference_wl);

cal_mua=water(index)*(1-alpha)+alpha*ua_hc(index); % this should equal to reference_mua

CSF_mua=(1-Arachnoid_fraction)*(water(:)*(1-alpha)+alpha*ua_hc)+Arachnoid_fraction*Arachnoid_mua;
CSF_mua=[wl CSF_mua];