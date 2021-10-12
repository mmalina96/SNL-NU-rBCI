clear all;
clc;
%% importing stuff, defining electrodes
cd C:\Users\micha\Documents\SNL_NU\FromMukta
load 'ChanLocs_128.mat'
nchan = 128;

%putting electrodes locs into matrices
allx = zeros(1,nchan);
ally = zeros(1,nchan);
allz = zeros(1,nchan);
for loc = 1:nchan
    x = locs(loc).X;
    y = locs(loc).Y;
    z = locs(loc).Z;
    if isempty(x)<1
        allx(1,loc) = x;
        ally(1,loc) = y;
        allz(1,loc) = z;
    end
end

central = [39 6 104 105 10 44 77 78 49 15 113 114 20 54 87 88 58 26 121 122];
central = sort(central);

periph = [71 46 22 124 32 128 96 125 28 117 17 41 66 65 2 36 7 45 16 55 ...
    27 64 30 31 29 60 23 51 12 42 3 33 1 108 81 91 123 95 127 126 92 118 ...
    82 72 99];
periph = sort(periph);


%% calculating euDist mat, nearest electrodes
electrodeLoc = [allx; ally; allz]';
distVec = zeros(128,128);
for j = 1:128
    for i = 1:128
        currentLoc = [allx(i); ally(i); allz(i)]';
        distVec(j,i) = sqrt((electrodeLoc(j,1) - currentLoc(1,1))^2 ...
            + (electrodeLoc(j,2) - currentLoc(1,2))^2 ...
            + (electrodeLoc(j,3) - currentLoc(1,3))^2);
    end
end
distVecNorm = normalize(distVec, 'range', [0,1]); %normalized moving fwd


% generating separate euDist maps for central and peripheral electrodes
for i = 1:length(central)
    electrode = central(i);
    centralLocs(i,:) = distVecNorm(electrode,:);
end

for i = 1:length(periph)
    electrode = periph(i);
    periphLocs(i,:) = distVecNorm(electrode,:);
end

% nearest electrode bit
[centralNeighborSort,CNindex] = sort(centralLocs, 2);
[periphNeighborSort, Pindex] = sort(periphLocs, 2);


%% generating Sim ICs
%cvErrMat = zeros(50,7,10);
% spatially relevant for sampling
periphICsTest = zeros(45,128,10);
centralICsTest = zeros(20,128,10);

% central case
centralICs = zeros(20,128);
keepMat = zeros(9000,128);
d = [1:20:9001]; %this will determine size of training mat
dist=0.1:0.1:0.3;
thresh = [0.1:0.2:0.9];
k = 0:0.2:1;
val = [2:2:8];

iterations=1;
cvErrMat_mix=zeros(1,size(dist,2)*size(k,2)*size(thresh,2)*6);
cvErrMat_pure=zeros(size(dist,2),size(k,2));
% =1;
wpcount=1;

w=1;
% cvdist
% cvk
% cvthresh
% cvmix
w_mix=zeros(size(dist,2)*size(k,2)*size(thresh,2)*length(val),128);
w_mixNorm=zeros(size(dist,2)*size(k,2)*size(thresh,2)*length(val),128);
w_mixACR=zeros(size(dist,2)*size(k,2)*size(thresh,2)*length(val),128);
w_mixEMG=zeros(size(dist,2)*size(k,2)*size(thresh,2)*length(val),128);
w_pure=zeros(size(dist,2)*size(k,2),128);
% for m=1:iterations
index_mix=zeros(size(dist,2)*size(k,2)*size(thresh,2)*length(val),4);
for dinc=1:length(dist)
    for iter_mat = 1:size(d,2)-1
        for i = 1:20
            for j = 1:128
                if centralLocs(i,j) == 0
                    centralICs(i,j) = 1;
                elseif  centralLocs(i,j) < dist(dinc) % what if you came up with a vector
                    % that was like (0.1:0.1:0.5) to rand determine the how far
                    % out the effected electrodes go. Might make generalization
                    % better
                    centralICs(i,j) = exp(-(centralLocs(i,j)).^2); % normal, test err = 1.1 %
                else
                    centralICs(i,j) = 0;
                end
            end
        end
        keepMat(d(iter_mat):(d(iter_mat+1)-1),:) = centralICs;
    end
    
    % peripheral case
    periphICs = zeros(45,128);
    throwMat = zeros(9000,128);
    e = [1:45:9001]; % this will determine size of training mat
    for iter_mat = 1:size(e,2)-1
        for i = 1:45
            for j = 1:128
                if periphLocs(i,j) == 0
                    periphICs(i,j) = 1;
                elseif  periphLocs(i,j) < 2*dist(dinc)
                    periphICs(i,j) = exp(-(periphLocs(i,j)).^2); % normal
                else
                    periphICs(i,j) = 0;
                end
            end
        end
        throwMat(e(iter_mat):(e(iter_mat+1)-1),:) = periphICs;
    end
    pureICs = [repmat(centralICs,10,1), ones(1,200)'; repmat(periphICs,10,1), (-1.*ones(1,450)')];
    
    % to properly model the negatives and positives, should I make whole rows
    % either all negative or all positive, or whole mats, or just individual
    % entries or what?
    
    
    %% pure model
    % general model: x(hEEG) + (1-x)EMG + k*rand(1,128)
    % x, k testable variable
    
    % pick random rows to choose from
    
    keepLabel = ones(9000,1);
    throwLabel = -1 .* ones(9000,1);
    dataMat_pure = [keepMat keepLabel; throwMat throwLabel];
    %dataMat = [dataMat; mixIC];
    % dataMat_mix = [mixIC];
    
    %add noise term
    % noise = 0.005 .* rand(18000,128);
    % noise = [noise zeros(18000,1)];
    % dataMat = dataMat + noise;
    
    % k=0.1;
    nrows=size(dataMat_pure,1);
    for i = 1:length(k)
        
        noise_pure = k(i) .* (rand(nrows,128)-0.5);
        noise_pure = [noise_pure zeros(nrows,1)];
        dataMat_pure = dataMat_pure + noise_pure;
        
        
        change_ind_p=randperm(nrows);
        dataMat_chk_p=dataMat_pure(change_ind_p,:);
        train_data_p=dataMat_chk_p;
        test_data=testDataLoad; % from characterize ICs script
        
        feat_train_p=train_data_p(:,1:128);
        y_train_p=train_data_p(:,129);
        feat_test_p=test_data(:,1:128);
        y_test=test_data(:,129);
        
        w_p=linclf(feat_train_p,y_train_p);
        
        y_pred_train_p=feat_train_p*w_p;
        y_pred_test_p=feat_test_p*w_p;
        
        y_pred_train_p(y_pred_train_p<0)=-1;
        y_pred_train_p(y_pred_train_p>=0)=1;
        
        y_pred_test_p(y_pred_test_p<0)=-1;
        y_pred_test_p(y_pred_test_p>=0)=1;
        
        %     fprintf("Training err:%f",abs(sum(y_pred_train_p-y_train_p))*100/size(y_pred_train_p,1));
        %     fprintf("\nTesting err:%f",abs(sum(y_pred_test_p-y_test))*100/size(y_pred_test_p,1));
        
        [tst_cverr_p,tr_cverr_p,w_pure(wpcount,:)]=crossval_linclf(dataMat_chk_p,10);
        
        cvErrMat_pure(dinc,i)=mean(tst_cverr_p);
        wpcount=wpcount+1;
        
    end
    
    
end
%% mix begins here
cd C:\Users\micha\Documents\SNL_NU\FromMukta\intentEMGs
load('acrMixMat');
load('emgMixMat');
cd C:\Users\micha\Documents\SNL_NU\FromMukta
wmcount=1;
%mixICs = zeros(9800,129,4);
for dinc=1:length(dist)
    for i=1:length(k)
        for t = 1:length(thresh)
            for j = 1:length(val)
                for thing = 1:9000
                    x = rand;
                    hEEGrow = randperm(20,val(j));
                    EMGrow = randperm(45,val(j));
                    mixICs(thing,1:128,j) = sum(((x.*(centralICs(hEEGrow,:)))./val(j)) + (((1-x).*(periphICs(EMGrow, :)))./val(j))); %+ k.*rand(1,128)];
                    if x <= thresh(t)
                        mixICs(thing,129,j) = -1;
                    else
                        mixICs(thing,129,j) = 1;
                    end
                end
            end
            ICs = randperm(9000,3000);
            mixICs2 = [mixICs(ICs,:,1)];
            mixICs4 = [mixICs(ICs,:,2)];
            mixICs6 = [mixICs(ICs,:,3)];
            mixICsCong = [mixICs2; mixICs4; mixICs6];
            mixICs(1:9000,:,4) = mixICsCong;
            mixICs(9001:9650,:,1) = pureICs;
            
            
            % adding intent EMG information
            mixICsACR = mixICs;
           
            
            mixICsACR(9651:10418,1:129,1) = acrMixMat;
            mixICsACR(9651:10418,1:129,2) = acrMixMat;
            mixICsACR(9651:10418,1:129,3) = acrMixMat;
            mixICsACR(9651:10418,1:129,4) = acrMixMat;
            
            mixICsEMG = mixICs;
            mixICsEMG(9651:10418,1:129,1) = emgMixMat;
            mixICsEMG(9651:10418,1:129,2) = emgMixMat;
            mixICsEMG(9651:10418,1:129,3) = emgMixMat;
            mixICsEMG(9651:10418,1:129,4) = emgMixMat;
            
            
            count=1;
            %             cvErrMat(dinc,i,y,count) = dist(dinc);
            
            for var=1:length(val)
                index_mix(w,:)=[dist(dinc),k(i),thresh(t),val(var)];
                dataMat_mixNorm = mixICs(:,:,var);
                dataMat_mixACR = mixICsACR(:,:,var);
                dataMat_mixEMG = mixICsEMG(:,:,var);
               [errNorm,w_mixNorm(wmcount,:)]=mixGraph(dataMat_mixNorm,k(i));
                [errACR,w_mixACR(wmcount,:)]=mixGraph(dataMat_mixACR,k(i));
                [errEMG,w_mixEMG(wmcount,:)]=mixGraph(dataMat_mixEMG,k(i));
                
                cvErrMat_mixNorm(1,wmcount)=errNorm;
                cvErrMat_mixACR(1,wmcount)=errACR;
                cvErrMat_mixEMG(1,wmcount)=errEMG;
                
                %             fprintf("%d\n",val(var));
                count=count+1;
                wmcount=wmcount+1

            end
            
        end
    end
end

%% create labelMat
labelMat = classifyICs(w_mix);

% creating new mat to determine rows (parameter spaces) that share label
% outputs
remove_cols=cell(size(labelMat));
for lcount=1: size(labelMat,1)
    diffMat = labelMat{lcount,1};
    diffMat(diffMat == -1) = 0;
    
    %unique func should return the unique rows
    [C, indexC] = unique(diffMat, 'rows');
    
    %this will tell us all combinations of rows
    diffVec = nchoosek(1:size(w_mix,1),2);
    
    for i = 1:length(diffVec)
        diffVal(i,1) = sum(abs(diffMat(:,diffVec(i,1))-diffMat(:,diffVec(i,2))));
        diffVal(i,2)= diffVec(i,1);
        diffVal(i,3) = diffVec(i,2);
    end
    %     diffMat(end,end)
    diffVal = sortrows(diffVal, 1);
    rem=find(diffVal(:,1)==0); %first column store difference between classifiers
    remove_cols{lcount,1}=zeros(size(rem,1),3);
    remove_cols{lcount}(:,1)=rem;
    remove_cols{lcount}(:,2:3)=diffVal(rem,2:3);
end

%% least sqaure linear classifiers

% for m=1:10
% % y = [0:0.1:0.9];
%     keepLabel = ones(9000,1);
%     throwLabel = -1 .* ones(9000,1);
%     dataMat_pure = [keepMat keepLabel; throwMat throwLabel];
%     %dataMat = [dataMat; mixIC];
%     % dataMat_mix = [mixIC];
%
%     %add noise term
%     % noise = 0.005 .* rand(18000,128);
%     % noise = [noise zeros(18000,1)];
%     % dataMat = dataMat + noise;
%
%     k = 0.1;
%     % k=0.1;
%     nrows=size(dataMat_pure,1);
%     for i = 1:length(k)
%
%         noise_pure = k .* (rand(nrows,128)-0.5);
%         noise_pure = [noise_pure zeros(nrows,1)];
%         dataMat_pure = dataMat_pure + noise_pure;
%
%
%         change_ind_p=randperm(nrows);
%         dataMat_chk_p=dataMat_pure(change_ind_p,:);
%         train_data_p=dataMat_chk_p;
%         test_data=testDataLoad; % from characterize ICs script
%
%         feat_train_p=train_data_p(:,1:128);
%         y_train_p=train_data_p(:,129);
%         feat_test_p=test_data(:,1:128);
%         y_test=test_data(:,129);
%
%         w_p=linclf(feat_train_p,y_train_p);
%
%         y_pred_train_p=feat_train_p*w_p;
%         y_pred_test_p=feat_test_p*w_p;
%
%         y_pred_train_p(y_pred_train_p<0)=-1;
%         y_pred_train_p(y_pred_train_p>=0)=1;
%
%         y_pred_test_p(y_pred_test_p<0)=-1;
%         y_pred_test_p(y_pred_test_p>=0)=1;
%
%         %     fprintf("Training err:%f",abs(sum(y_pred_train_p-y_train_p))*100/size(y_pred_train_p,1));
%         %     fprintf("\nTesting err:%f",abs(sum(y_pred_test_p-y_test))*100/size(y_pred_test_p,1));
%
%         [tst_cverr_p,tr_cverr_p]=crossval_linclf(dataMat_chk_p,10);
%         cvErrMat(:,1,m) =k;
%         cvErrMat(:,2,m) = mean(tst_cverr_p);
%     end
%
% %     count=3;
% %     for var=1:length(val)
% %         dataMat_mix = mixICs(:,:,var);
% %         cvErrMat(:,count,m)=mixGraph(dataMat_mix);
% %         fprintf("%d\n",val(var));
% %         count=count+1;
% %     end
% end

% cvErrMatMean = mean(cvErrMat,3);
%
% figure();
% plot(cvErrMatMean(:,1),cvErrMatMean(:,2));
% hold on;
% plot(cvErrMatMean(:,1),cvErrMatMean(:,3));
% hold on;
% plot(cvErrMatMean(:,1),cvErrMatMean(:,4));
% hold on;
% plot(cvErrMatMean(:,1),cvErrMatMean(:,5));
% hold on;
% plot(cvErrMatMean(:,1),cvErrMatMean(:,6));
% hold on;
% plot(cvErrMatMean(:,1),cvErrMatMean(:,7));
% hold on;
% plot(cvErrMatMean(:,1),cvErrMatMean(:,8));
% xlabel("Noise amplitude");
% ylabel("Average CV test error");
% legend("Pure+noise","Mixed+noise2","Mixed+noise4","Mixed+noise6","Mixed+noise8","Mixed+noise10", "Mixed+noise Combined")
% title("CV plots");

%