function [result,w_mix] = mixGraph(dataMat_mix,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% least sqaure linear classifiers

%dataMat = [dataMat; mixIC];
% dataMat_mix = [mixICs];

%add noise term
% noise = 0.005 .* rand(18000,128);
% noise = [noise zeros(18000,1)];
% dataMat = dataMat + noise;

% k = 0.1:0.1:1;
nrows_mix=size(dataMat_mix,1);
result=zeros(length(k),1);
% thresh = [0:0.1:0.9];
% for i = 1:length(k)
    
    
    %%%%%mix begins here
    
    noise_mix= k .* (rand(nrows_mix,128)-0.5);
    noise_mix = [noise_mix zeros(nrows_mix,1)];
    dataMat_mix = dataMat_mix + noise_mix;
    
    change_ind_m=randperm(nrows_mix);
    dataMat_chk_m=dataMat_mix(change_ind_m,:);
    train_data_m=dataMat_chk_m;
    %     test_data=testDataLoad; % from characterize ICs script
    
    feat_train_m=train_data_m(:,1:128);
    y_train_m=train_data_m(:,129);
    %     feat_test_m=test_data(:,1:128);
    %     y_test=test_data(:,129);
    
    w_m=linclf(feat_train_m,y_train_m);
    
    y_pred_train_m=feat_train_m*w_m;
    %     y_pred_test_m=feat_test_m*w_m;
    
    y_pred_train_m(y_pred_train_m<0)=-1;
    y_pred_train_m(y_pred_train_m>=0)=1;
    
    %     y_pred_test_m(y_pred_test_m<0)=-1;
    %     y_pred_test_m(y_pred_test_m>=0)=1;
    
    %     fprintf("Training err:%f",abs(sum(y_pred_train_m-y_train_m))*100/size(y_pred_train_m,1));
    %     fprintf("\nTesting err:%f",abs(sum(y_pred_test_m-y_test))*100/size(y_pred_test_m,1));
    
    
    [tst_cverr_m,tr_cverr_m,w_mix]=crossval_linclf(dataMat_chk_m,10);
    
    %     cvErrMat(i,1) = k(i);
    result = mean(tst_cverr_m);
end


