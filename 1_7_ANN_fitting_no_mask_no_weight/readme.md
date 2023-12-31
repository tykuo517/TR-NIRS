# ANN fitting
Generate target spectra, or use the measured target spectra, and fitting them.

---
## Prepare
* An ANN for each subject that inputs the OPs of each layers and output the reflectance spectra for each SDSs.  
* The boundary of each tissue parameter.  
* The calibrated target spectra.  

---

## Background

To make the optical parameters (OPs) meaningful continuous in each wavelength, I use tissue parameters to describe the OPs of each layer.  
The tissue parameters including:
* For mua: the concentration of each absorber, e.g., hemoglobin (hemoglobin concentration (hc) and tissue oxygen saturation (StO2) are used to calculate the concentation of oxyhemoglobin and deoxyhemoglobin) and melanin.
    * The true mua is calculated by multiply the absorption coefficient of each absorber by the their concentraiton.
* For mus: A and K are used to describe the scattering coefficient.
    * $\mu_s'(\lambda)=A\times\lambda^{-k}$
    * $\mu_s=\mu_s'/(1-g)$

---

## 1. Make initial value database
Because the fitting needs a initial value to begin with, we need to use a database to choose which initial value(s) is better. The database(DB) contains many sets of tissue parameters (which describe the OPs for each layer) and the corresponding spectra. Before fitting, the program will compare all spectra in the database to the target spectrum, and choose the tissue parameters with the smallest spectral error as the initial value of fitting.

### 1.1. Find AK boundary

Before make the DB, since the mus is describe by A and K, which is a nonlinear relationship, we should make sure the calculated mus is not exceed the mus boundary(of the ANN).


Use `S0_find_AK_range.m` to find the boundary.

You should give the name of the ANN models and the dir of the ANN models, also the g value and wavelength range you will used to fitting.
```matlab=11
subject_name_arr={'ZJ','WW2','YF','YH','WH2','KB','SJ','BT','SC'}; % the name of the subjects
g=0.9; % g for L1,2,4
target_spec_wl=(700:900)';

model_folder='model_arrange'; % the folder of the models
```
Then the program will load the ANN model, find the mus boundary of the ANN, and auto calculat the range of A and K.

### 1.2. Make database

Use `S1_make_database.m` to make the database.

You should give the dir of ANN, the name of ANN models.  
`fitting_wl` is the wavelength to generate the spectra corresponding to the initial values. Since the wavelength we used to fitting may be different from time to time, the spectra in the DB should be fine(dense in wavelength) enough to let the interpolation work. Also, we should also check the wavelength we will used to fitting in `fitting_wl_file`, to check if the randomly generated mua will exceed the ANN boundary in the fitting wavelength.  
You should also give the `mua_param_Lbound` and `mua_param_Ubound` of each tissue parameters.

```matlab=13
model_dir='model_arrange'; % the folder of the models
subject_name_arr={'ZJ','WW2','YF','YH','WH2','KB','SJ','BT','SC'}; % the name of the subjects
num_SDS=7; % the number of SDS, in the ANN

num_init_spec=50000; % number of random generated spectrum
fitting_wl=(700:4:900)'; % the wavelenth to generate the spectrum

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

output_dir='initValue_database_2'; % the folder to save the output database

% the range for the parameters
%                      hc_1    sto2_1  hc_2    sto2_2 hc_4     sto2_4  mel
mua_param_Lbound=     [1       0.3     1       0.3    50       0.3     0];
mua_param_Ubound=     [150     1       150     1      200      1       0.003 ];
```

The generated DB will be one file for each subject.  
![](https://i.imgur.com/2f9Z8T8.png)

### 1.3. Check the DB

If you had changed the fitting wl, you should use `S1_2_database_refine.m` and `S1_3_check_refined_database.m` to make sure the DB will not exceed the ANN boundary. If there are initial value that exceed the boundary, the program will auto replace them to fix the problem.

## 2. Generate testing spectra

To generat the testing spectrum, you should use `T2_generate_target_paramAnswer.m`.  

In the program, you should give:
* `num_anser_to_generate`: the number of random parameter to generate. For example, there will be 15 random answer to fitting for each subejct.  
* `num_error_to_generate`: the number of error to add to each testing spectrum. Since the invivo measurement will have some noise, we should add these noise to the target spectrum(otherwise will have no noise, since it's generated by ANN). The first one will be no error. For example, `num_error_to_generate=15` means there will be 1 spectrum with no noise and 14 spectra with error.
* `SDS_error_arr` is the CV of each SDS when you measure the phantom or subject.
* `use_slope_mode` is used if you want to generate the noise representing the invivo noise. Since I found that the invivo nosie will be differnet in each wavelength, the flag `use_slope_mode=1` will let the program generate different error for each wavelength. Otherwile, `use_slope_mode=0` will generate the same noise for each wavelength to represent the noise while measure phantom, which is the system noise.
* `small_noise_ratio` will generate additional small noise in each wavelength.

```matlab=14
subject_name_arr={'ZJ','WW2','WH2','YF','KB','YH','SJ','BT','SC'};
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
% SDS_error_arr=[2 2.3 0.75 4 1 2]; % in %, the error (CV) of each SDS
SDS_error_arr=[3 4.2 5.1 5.2 5.4 12.1]; % in %, the error (CV) of each SDS
use_slope_mode=1; % if =1, use slope noise instead of constant noise
small_noise_ratio=0.05;
```

After generate the target spectra, you can use `T2_2_find_spec_addedNoise.m` to calculate the CV of the generated spectra to find if the CV is near the value you set. `T2_3_plot_generated_spec.m` can be used to plot the generated spectra(with differnet noise).
![](https://i.imgur.com/Q8CkKgn.jpg)

## 3. Fitting testing target spectra

### 3.1. Fitting
Use `T3_main_fitting_generated_target.m` to fitting the generated spectra.

You should give:
* `times_to_fitting`: how many initial values are used to fitting. Use more than only 1 initial value to prevent local minima.
* `input_dir`: the folder containing the testing spectra you just generated.
* `SDS_choosed`: Which SDS you choose to fitting. The program will only calculate the error of these SDSs and use the error to optimize the OPs.

```matlab=13
subject_name_arr={'ZJ','WW2','WH2','YF','KB','YH','SJ','BT','SC'};
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value
input_dir='test_fitting_2021-01-17-16-47-46'; % the folder containing the generated target spectrum

num_SDS=6; % how many SDS are in the target spectrum
SDS_choosed=[1 2 3 4 5 6]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength
```

After the fitting is done, you can use `T3_2_check_fitting_complete.m` to chech of all the fitting are done. Sometimes there will be some error for some initial value, make sure you had run `S1_2_database_refine.m` if you had change the fitting wl. Otherwise, you can use `T3_main_fitting_generated_target.m` to fitting again, the program will skip the target which had been done automatically.

### 3.2. Arrange the fitting result

Use `T4_fitting_result_arrange.m` to arrange the fitting result.  

You should give which SDS used to fitting and the `fitting_dir` is the folder that `T3_main_fitting_generated_target.m` create automatically for each SDS combination.
```matlab=19
fitting_dir='fitting_SDS1234'; % the fitting folder

SDS_choosed=[1 2 3 4]; % the SDSs you used to fitting
```
The program will arrange the fitting result of multiple initial values, and find which fitting have the smallist error for choosed SDSs.

There will be one .csv file for each spectra
![](https://i.imgur.com/xVuzNzH.png)
Some figure are also plot, including the fitted spectra for each initial value and the fitted OPs.
![](https://i.imgur.com/LH1WPqq.png)

### 3.3. Result analysis

You can use `T6_plot_closeError_spec_OP.m` to plot some fitting result which have the smallist fitted spectral error.
![](https://i.imgur.com/4QFZfWI.png)

You can use `T8_cal_OP_spec_error_corr.m` to plot the correlation coefficient of the spectral error and OPs error.
![](https://i.imgur.com/w1ui5we.png)

You can use `T9_cal_bestFitting_opError_percentile.m` to calculate the OP error percentile for each rank's fitting result.
![](https://i.imgur.com/d29euxI.png)

Using the result above, I think that if I choose the fitting result with smallist spectral error, I can almost get the smallist OPs error. So you can use `T10_cal_mean_OPerror.m` to calculate the mean fitted OP error for all target spectra.

If you had fitting the invivo subject spectra, you can use `T12_save_CI_for_subject.m` and `T13_save_errorPDF_for_subject.m` to copy the informaiton of OP error confidence interval or the PDF of the error distribution. Since each subject may use different SDS set to fitting, the confidence interval or error PDF may be different.

## 4. Fitting invivo spectra

### 4.1. Fitting

You should first put the calibrated target spectra in one folder:
![](https://i.imgur.com/6cV1smK.png)  

Use `S2_main_fitting_spec.m` to fitting the invivo spectra.  

The parameter to give:
* `target_dir`: the folder containing the calibrated target spectra.
* `SDS_choosed`: Which SDS you choose to fitting. The program will only calculate the error of these SDSs and use the error to optimize the OPs.
*  `times_to_fitting`: how many initial values are used to fitting. Use more than only 1 initial value to prevent local minima.

```matlab=13
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
model_name_arr={'ZJ','WW2','WH2','YF','KB'}; % the name of ANN model corresponding to each target spec

output_folder='fitting_SDS346_3'; % the name of output main folder

num_SDS=6; % how many SDS are in the target spectrum
SDS_choosed=[3 4 6]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

times_to_fitting=20; % number of init value used to fitting

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

model_dir='model_arrange'; % the folder containing the arranged ANN file
target_dir='input_target_2'; % the folder containing the target spectrum
```

### 4.2. Arrange the fitting result

Use `S3_fitting_result_arrange.m` to arrange the fitting result.

### 4.3. Fitting result analysis

You can use `S4_compare_fittedParam_toError.m` to plot the scattering plot of the fitted tissue parameters and spectral error.
![](https://i.imgur.com/aItVkvT.png)

You can use `S5_plot_closeError_spec_OP.m` to plot some fitting result which have the smallist fitted spectral error.
![](https://i.imgur.com/vMAXKfS.png)

### 4.4. Choose the fittid OPs

To prevent multiple solution, you can use `S6_choose_save_fitted_OP.m` to choose which answer is the best answer, or there are multiple solution.

You should give the parameters:
* `error_range_threshold`: the threshold of error. If the difference between two fitting's spectral error is smaller than this threshold, then they are both consider possible solutions. You may use the average system noise as the threhsold.
* `toOutput_subject_index`: which subject to output. Because this program will choose between multiple solution, each time only process one subject.
* `fix_toOutput_rank`: If you already know which rank of the fitting result if the suitable solution(s), you can give this parameter and skip the choose peocess. Otherwise, the program will choose the suitable solution.

```matlab=13
error_range_threshold=0.02; % The error difference smaller than this value will be consider multiple solution
output_dir='fitted_result'; % the folder to output the fitted OP
to_output_wl=(700:890)'; % The wavelength to output the fitted OP

fitting_dir={'fitting_SDS1234_3','fitting_SDS2345_3','fitting_SDS123456_3','fitting_SDS12345_3','fitting_SDS346_3'}; % The fitted folder for different SDS combination
fitted_OP_error=[13.2 5.7 26.3 10.7 23.1 56.2 % the error of fitted OP for each SDS combinations
                   15.4	8.1 26.4 14   24.5 51.7
                   12.4	5.4 24.7 10   19.8 40
                   12.5 5.4 26.5 10.5 22.9 44.9
                   15.7 7.2 30.4 14.0 22.4 55.8];
compensate_SDS_index_arr={[5 6],[1 6],[],[6],[1 2 5]}; % the SDSs didn't used to fitting

sbj_fitting_dir_arr=[3 5 1 4 1]; % which fitting folder does this subject use
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
model_name_arr={'ZJ','WW2','WH2','YF','KB'}; % the name of ANN model corresponding to each target spec

toOutput_subject_index=1; % the index of the subject to output
% fix_toOutput_rank=[2 8]; % skip the chossing process, output these ranks
fix_toOutput_rank=[]; % skip the chossing process, output these ranks
```

The program will first compare the SDSs choosed to fitting, if the distance between the smallest error and other error is larger than the threshold, then the smallist error will be choosed as the only one solution.  
If there are other result that have small enough error, then we will compare the error of other SDSs(which not be used to fitting). As the fitting SDS sets, the program will find the smallest error, and find if there are other result which also have small error(difference smaller than threshold).  

The choose process will be plot.  
* Only one solution:
![](https://i.imgur.com/QAcpZa9.png)

* Two solution:
![](https://i.imgur.com/A9PQNaF.png)

The OP and fitted spectrum of the choosed result(s) will be in the output dir.
![](https://i.imgur.com/ShLksua.png)

You can use `S7_plot_fittedOP_compareToLiterature.m` to plot the choosed results and compare them to the literature value.

* `add_OP_error`: if =1, the error will be present in shaded errorbar.
* `use_CI_as_error`: if =1, the plot will use the confidence interval as the error. Otherwise, it will use the mean error.
* `CI_index`: plot the first confidence interval (68%) or the second confidence intervel (95%).
    * P.S. If you want to use confidence interval as the error, you should run `T12_save_CI_for_subject.m`.
```matlab=16
to_plot_layer=[1 1 0 1]; % if plot scalp, skull, CSF and GM
add_OP_error=1; % if =1, add the OP error to the fitted result
use_CI_as_error=1; % if =1, use confidence interval rather than CI to plot
CI_index=2; % if use_CI_as_error, the confidence interval to plot
```

The plot result:
![](https://i.imgur.com/ooIAFad.png)

You can also compare individual subject's result to the literature value. This may help you to choose if all the multple solution is reasonable.
![](https://i.imgur.com/TradrEA.png)

If you find some of the multiple solution is more reasonable than other solution, you can use `S7_3_process_multiSol.m` to process the original choosed data, let the new choosed data as the result.

* `to_process_sbj_index`: which subject to process
* `to_save_multiSol_index`: the index of the solution you consider is more suitable. If the second originally choosed result is better, you should type 2 instead of the true rank if it.
```matlab=14
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
to_process_sbj_index=3; % the index of the subject to process
to_save_multiSol_index=[1]; % the array containing the index of the solution(s) to save
```

And you can use `S7_4_plot_multiSol_processed_spec.m` to replot the choosed spectra for the multisolutin subject.
![](https://i.imgur.com/TMU2bug.png)

## 5. Output the fitting OPs

### 5.1. Output OPs

You can use `S8_output_fitted_OP.m` to simply output the choosed fitting result.

You can also use `S8_2_output_fitted_OP_random.m` to add noise (fitted OP error) to the OPs while output the fitted OPs.  
This program will use Gaussian distribution to sample the nosie with mean value = mean fitted OP error (found using testing target spectra).

* `num_random_OP`: how many times of random noise to generate.
* `to_error_OP`: which OP needs to add error
``` matlab=17
num_random_OP=15; % number of random nosie to added to the OP

add_OP_error=1; % if =1, add the OP error to the fitted result
to_error_OP=[1 2 3 4 7 8]; % the OPs needs to add error, 1 = mua1, 2 = mus1, 3 = mua2 ......
```

You can also use `S8_2_output_fitted_OP_random_PDF.m` to add nosie to the OPs.  
This time, the noise is generated according to the error PDF (found using testing target spectra)  
P.S. You should run `T13_save_errorPDF_for_subject.m` first to save the error PDF for each subject.  

* `num_random_OP`: how many times of random noise to generate.
* `to_error_OP`: which OP needs to add error
```matlab=17
num_random_OP=30; % number of random nosie to added to the OP

to_error_OP=[1 2 3 4 7 8]; % the OPs needs to add error, 1 = mua1, 2 = mus1, 3 = mua2 ......
```

### 5.2. Compare ANN spec to MC spec

After you use the fitted OP to run MC simulation, you can use `S9_compare_fitted_spec_ANN_MC.m` to compare the fitted spectrum generated by ANN and MC simulation.
![](https://i.imgur.com/N5WXPH0.png)
