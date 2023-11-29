# MCX ANN training
Use the pre-generated training data to train the ANN

---

## Prepare
* The pre-simulated training data containing the mua, mus for each layer, and the reflectance for each SDS.

---

## Simulation steps

1. In `S1_main_training_allSDS.m`, set the model name also the path to the training data.  
    ```matlab=12
    subject_name_arr={'ZJ','WW','YF','YH','WH','KB','BT','SJ','SC'}; % the name of the subject
    input_mother_dir='../20200622_MCX_lkt_trainingDataGen'; % the folder of the lookup table forward training data
    input_dir_arr={'ZJ_2020-12-21-23-54-25','WW_2021-01-04-07-36-02','YF_2020-12-22-00-05-23','YH_2020-12-22-01-47-27','WH_2021-01-05-07-15-12','KB_2020-12-22-00-41-56','BT_2020-12-22-01-35-21','SJ_2020-12-22-08-54-09','SC_2020-12-22-01-46-02'}; % the lut folder of each subject
    ```
    
    
2. Set the ANN setting:  
    ```matlab=22
    num_SDS=7;
    do_normalize=1; % if =1, do normalization to spec and param
    max_sequence_length=50000; % how many data sets in one batch, to prevent use too much GPU memory

    testing_rate=0.15; % the ratio of testing data
    validation_rate=0.1; % the ratio of validation data
    ```
    * num_SDS: how many SDS to predict.  
    * max_sequence_length: how many sets of mua/mus/reflectance to be in one group, should be smaller if your GPU have little RAM space.  
    
    The ANN training options:
    ```matlab=90
    options = trainingOptions('adam', ...
        'Shuffle','every-epoch', ...
        'ExecutionEnvironment','gpu', ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.95, ...
        'LearnRateDropPeriod',10, ...
        'ValidationData',validation_data, ...
        'MaxEpochs',1000, ...
        'ValidationFrequency',num_training_group, ...
        'ValidationPatience',20, ...
        'MiniBatchSize',128, ...
        'SequenceLength',max_sequence_length, ...
        'VerboseFrequency',1, ...
        'Plots','training-progress');
    ```
    
    The structure of ANN:
    ```matlab=105
    layers=[sequenceInputLayer(8,'Name','')
    %         fullyConnectedLayer(650)
        fullyConnectedLayer(850)
        leakyReluLayer
    %     reluLayer
    %         fullyConnectedLayer(450)
        fullyConnectedLayer(550)
        leakyReluLayer
    %     reluLayer
        fullyConnectedLayer(300)
        leakyReluLayer
    %     reluLayer
        fullyConnectedLayer(150)
        leakyReluLayer
    %     reluLayer
        fullyConnectedLayer(7)
        regressionLayer];
    ```
    
3. Run `S1_main_training_allSDS.m` to train the ANN. Matlab will show the following window automatically.  
    ![](https://i.imgur.com/1FYqxWK.png)  
    If you think the loss is not decay anymore, then you can stop the trainining manually by clicking the stop button.  
    ![](https://i.imgur.com/iCKzikM.png)  
    
4. The training result will be in a [name + date] folder, containing the trained ANN, the error record, also the backup file of `S1_main_training_allSDS.m` for easier tract the setting.  
    ![](https://i.imgur.com/KQcaDhh.png)  
    Also a error histogram show the goodness of the training.  
    ![](https://i.imgur.com/RvumegZ.png)  

5. Analysis the performance of the ANN.  

    You can use `S2_1_main_ANN_forward.m` to use the ANN to generate some forward spectrum, and compare it to the simulation.

    You can also use `S2_2_ANN_error_analysis.m` to analysis the where does the error come from  
    ![](https://i.imgur.com/R8KRE7N.png)

    Use `S2_3_rePlot_error_hist.m` to plor the error histogram (better looking one)  
    ![](https://i.imgur.com/5OGNdF6.png)

    Use `S2_4_plot_training_progress.m` to plot the training progress.  
    ![](https://i.imgur.com/xglIAyk.png)
    
6. Use `S3_save_model.m` to save the model into a `.mat` file.

    You can use a file to store the subject name and the corresponding ANN training folder.  
    ![](https://i.imgur.com/qvrZbdk.png)

7. You can use `S4_train_error_arrange.m` to print the error for each SDS and subject together.  
    ![](https://i.imgur.com/0YxI3Px.png)


