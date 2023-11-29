%{
Convert the tissue parameter into the optical parameters

Inputs:
param_arr: [A1 K1 A2 K2 A3 K3 hc_1 sto2_1 hc_2 sto2_2 hc_4 sto2_4]
do_white: if =1, set the mua to 0

Outputs:
mu: [mua1 mus1 mua2 mus2 mua3 mus3 mua4 mus4 mua5 mus5]
mua_arr: [mua1 mua2 mua3 mua4 mua5]

Benjamin Kao
Last update: 2020/12/25
%}

function [mu,mua_arr]=fun_param_to_mu(param_arr,do_white)
    global lambda;

    %% param
    plot_fig=0;
    
    %% main
    A_arr=param_arr(1:2:5);
    K_arr=param_arr(2:2:6);
    abs_arr=param_arr(7:13);
    
    param = containers.Map;

    % sclap
    param('mel_1') = abs_arr(7); % 0.000664;
    param('A_1') = A_arr(1);
    param('K_1') = K_arr(1);
    param('hc_1') = abs_arr(1);
    param('sto2_1') = abs_arr(2);
    param('g_1') = 0.9;
    
    % skull, from tin CS result
    param('A_2') = A_arr(2);
    param('K_2') = K_arr(2);
    param('hc_2') = abs_arr(3);
    param('sto2_2') = abs_arr(4);
    param('g_2') = 0.9;
    
    % CSF
    param('g_3') = 0.9;
    
    % gray matter, from tin CS result
    param('A_4') = A_arr(3);
    param('K_4') = K_arr(3);
    param('hc_4') = abs_arr(5);
    param('sto2_4') = abs_arr(6);
    param('g_4') = 0.9;

    mu=param_to_mu_spec(param);
    mua_arr=mu(:,1:2:10);
    
    if do_white==1
        mu(:,1:2:10)=0;
    end
    

    %% plot mu figure
    if plot_fig==1
        close all;
        figure;
        hold;
        for i=1:5
            plot(lambda,mu(:,2*i-1));
        end
        legend({'scalp' 'skull' 'CSF' 'gray matter' 'white matter'},'Location','best');
        title('\mu_a');

        figure;
        hold;
        for i=1:5
            plot(lambda,mu(:,2*i));
        end
        legend({'scalp' 'skull' 'CSF' 'gray matter' 'white matter'},'Location','best');
        title('\mu_s');
    end
end


function output = param_to_mu_spec(param)
    global lambda;
    
    %% parameter
    skull_add_offset=0;
    
    global Mu_coeff;

    %% calculate mua, mus for tissue layers
    
    % scalp
    mua_hc = param('sto2_1') .* Mu_coeff('hbo') + ( 1 - param('sto2_1') ) .* Mu_coeff('hb'); % the hc_1 is uM, so directly calculate the Mu_a
    mu_a1 = param('hc_1') * mua_hc + 0.02 * Mu_coeff('collagen') + 0.75 * Mu_coeff('water') + param('mel_1') * Mu_coeff('melanin');
    mu_s1 = param('A_1') * ( lambda / 810 ) .^ -param('K_1') ./ ( 1 - param('g_1') );

    % skull
    if skull_add_offset==1
        skull_offset=0.2218; % 1/cm
    else
        skull_offset=0;
    end
    mua_hc = param('sto2_2') .* Mu_coeff('hbo') + ( 1 - param('sto2_2') ) .* Mu_coeff('hb'); % the hc_2 is uM, so directly calculate the Mu_a
    mu_a2 = param('hc_2') * mua_hc + 0.02 * Mu_coeff('collagen') + 0.75 * Mu_coeff('water') + skull_offset;
    mu_s2 = param('A_2') * ( lambda / 810 ) .^ -param('K_2') ./ ( 1 - param('g_2') );

    % CSF
    mu_a3 = Mu_coeff('CSF_mua');
    mu_s3 = Mu_coeff('CSF_mus');
    
    % gray matter
    mua_hc = param('sto2_4') .* Mu_coeff('hbo') + ( 1 - param('sto2_4') ) .* Mu_coeff('hb'); % the hc_5 is uM, so directly calculate the Mu_a
    mu_a4 = param('hc_4') * mua_hc + 0.75 * Mu_coeff('water');
    mu_s4 = param('A_4') * ( lambda / 810 ) .^ -param('K_4') ./ ( 1 - param('g_4') );
    
    % white matter
    mua_hc = param('sto2_4') .* Mu_coeff('hbo') + ( 1 - param('sto2_4') ) .* Mu_coeff('hb'); % the hc_5 is uM, so directly calculate the Mu_a
    mu_a5 = param('hc_4') / 2 * mua_hc + 0.75 * Mu_coeff('water'); % the hc of white matter is about half of hc of gray matter
    mu_s5 = param('A_4') * 3 * ( lambda / 810 ) .^ -param('K_4') ./ ( 1 - param('g_4') ); % the mus of white matter is about 3 times of gray matter
        
    output = [mu_a1 mu_s1 mu_a2 mu_s2 mu_a3 mu_s3 mu_a4 mu_s4 mu_a5 mu_s5];
end