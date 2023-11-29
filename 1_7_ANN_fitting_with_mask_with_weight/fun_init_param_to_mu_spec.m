%{
Preload the mu specs to avoid repeat computation

Benjamin Kao
Last update: 2020/09/19
%}

function fun_init_param_to_mu_spec()
global Mu_coeff;
global lambda;

epsilon_folder = 'epsilon';
%% load absorption coefficient
epsilon_arr = {'collagen','melanin','hb','hbo','water','CSF_mua','CSF_mus'};
epsilon_file_arr = {'collagen_1_cm.txt' 'mel_1_cm.txt' 'Hb_1_cm_um_mua.txt' 'HbO_1_cm_um_mua.txt' 'water_1_cm.txt' 'CSF_mua_1_cm.txt' 'CSF_mus_1_cm.txt'};
Mu_coeff = containers.Map;
for i = 1 : length(epsilon_arr)
    temp_mu_spec = load(fullfile(epsilon_folder, epsilon_file_arr{i}));
    Mu_coeff(epsilon_arr{i}) = interp1(temp_mu_spec(:,1), temp_mu_spec(:,2), lambda); % interpolate to used wavelength
end
end