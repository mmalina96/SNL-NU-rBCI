tic
% signal1 = load('GE_handS004R04_mvcut.mat')
% signal = signal1.signal(:,1:128);
% force = signal1.signal(:,129:130);
% states=signal1.states;
% dirtyParams = dirty.parameters;;


% function [cleanDataNorm, cleanDataACR, cleanDataEMG] = cleanICs(dirtyData, w_mixNorm,w_mixACR, w_mixEMG)
% load('w_mixACR')
load('w_mixNorm')
% load('w_mixEMG')
trial=["04","05","06","07","08","09","10"];
%dirtyData has size (588600,128)

for outer = 1:length(trial)
    wmatfile=sprintf("icaw_ge4%s_65_115.mat",trial(outer));
    wMat = load(wmatfile);
    wMat = wMat.ica_weights;
    wMatInv = inv(wMat);
    
%     ACRclass = w_mixACR * wMat ; ACRclass(ACRclass >= 0) = 1; ACRclass(ACRclass < 0) = -1;
    normClass = w_mixNorm * wMat; normClass(normClass >= 0) = 1; normClass(normClass < 0) = -1;
%     EMGclass = w_mixEMG * wMat; EMGclass(EMGclass >= 0) = 1; EMGclass(EMGclass < 0 ) = -1;
    
    wMat_buffbuffNorm = zeros(128,128,360);
    wMat_buffbuffACR = zeros(128,128,360);
    wMat_buffbuffEMG = zeros(128,128,360);
    
    for var_acr = 1:size(normClass,1)
        
%         rowACR=ACRclass(var_acr,:);
%         st_indACR=find(rowACR==1);
%         wMat_buffACR=wMat;
%         wMat_buffACR(st_indACR,:)=0;
%         wMat_buffbuffACR(:,:,var_acr) = wMat_buffACR;
        
        rowNorm = normClass(var_acr,:);
        st_indNorm = find(rowNorm == 1);
        wMat_buffNorm = wMat;
        wMat_buffNorm(st_indNorm,:)=0;
        wMat_buffbuffNorm(:,:,var_acr) = wMat_buffNorm;
        
%         rowEMG = EMGclass(var_acr,:);
%         st_indEMG = find(rowEMG == 1);
%         wMat_buffEMG = wMat;
%         wMat_buffEMG(st_indEMG,:)=0;
%         wMat_buffbuffEMG(:,:,var_acr) = wMat_buffEMG;
        
    end
    
    
    hemi = {'Fz','F4','F6','FFC4h','FFC6h','FC2','FC4','FC6','FCC4h','FCC6h','C4','C6','CCP4h','CCP6h','CP4','CP6','force','force2'};
    homolog = {'F1', 'F3', 'F5', 'FFC3h', 'FFC5h', 'FC1', 'FC3', 'FC5', 'FCC3h', 'FCC5h', 'C3', 'C5', 'CCP3h', 'CCP5h', 'CP3', 'CP5', 'force', 'force2'};
%     channels={hemi,homolog}; %index1 is hemi, index2 is homo
    forchan=load('parameters');
    channels=forchan.parameters.ChannelNames.Value;
    homo_ind=find(ismember(channels,homolog));
    hemi_ind=find(ismember(channels,hemi));
    chanind=[hemi_ind homo_ind];
    
    for ch = 1: size(chanind,2)
        if ch==1
            disp("Beginning hemi")
        else
            disp("Beginning homoLog")
        end
        file_name = sprintf("GE_handS004R%s_mvcut.mat",trial(outer));
        signal1 = load(file_name);
        signal = signal1.signal(:,1:128);
        force = signal1.signal(:,129:130);
        states = signal1.states;
%         chan = channels{ch};
        
%         
%         signal1Norm.parameters.ChannelNames.Value=chan; %change depending on homoLog or homoLog
%         signal1ACR.parameters.ChannelNames.Value=chan; %change depending on homoLog or homoLog
%         signal1EMG.parameters.ChannelNames.Value=chan; %change depending on homoLog or homoLog
%         
        
        vafMatNorm = zeros(360,2);
        vafBadLogicalNorm = zeros(360,128);
        
%         vafMatACR = zeros(360,2);
%         vafBadLogicalACR = zeros(360,16);
%         
%         vafMatEMG = zeros(360,2);
%         vafBadLogicalEMG = zeros(360,16);
        
        
        signal1Norm.parameters = signal1.parameters;
        signal1Norm.parameters.ChannelNames.Value = signal1.parameters.ChannelNames.Value; %change depending on homoLog or homoLog
        signal1Norm.states = signal1.states;
        
%         signal1ACR.parameters = signal1.parameters;
%         signal1ACR.parameters.ChannelNames.Value=chan; %change depending on homoLog or homoLog
%         signal1ACR.states = signal1.states;
%         
%         signal1EMG.parameters = signal1.parameters;
%         signal1EMG.parameters.ChannelNames.Value=chan; %change depending on homoLog or homoLog
%         signal1EMG.states = signal1.states;
        
%         signal1.parameters.ChannelNames.Value=chan;
        
        %%
        for var_normc = 1:size(normClass,1)
            
            %norm
            cleanDataNorm = signal' - (wMat_buffbuffNorm(:,:,var_normc) * wMat * signal');
%             cleanDataNorm = cleanDataNorm(chanind(1:end-2,ch),:)';
            signal1Norm.signal=[cleanDataNorm', force];
        if ch==1
            %pathname = fileparts('Results/');
            pathname = 'C:\Users\micha\Documents\SNL_NU\fromPrashanth\EMG_artifact_detect_update\EMG artifact detection\Results'; %changed this last night for it to run, local path
            [vafNorm, badLogicalNorm, bd_hemi] = cleanVAF(signal1Norm,ch,chanind(:,ch));
            savebadchans = append(file_name,"bad_chans.mat");
            matfile_badchans = fullfile(pathname,savebadchans);
            save(matfile_badchans,'bd_hemi');
        else
%             file_hemi = sprintf("GE_handS004R%s_1",trial(outer)); %file_count isn't defined anywhere
            name_hemi = append(file_name,"bad_chans.mat");
            ld_hemi = load(name_hemi);
            bd_hemi = ld_hemi.bd_hemi;
            [vafNorm,badLogicalNorm] = cleanVAF(signal1Norm,ch,chanind(:,ch),bd_hemi);
        end
            vafMatNorm(var_normc,1) = mean(vafNorm);
            vafMatNorm(var_normc,2) = sum(badLogicalNorm);
            vafBadLogicalNorm(var_normc,:) = badLogicalNorm;
            
            %acr
%             cleanDataACR = signal' - (wMat_buffbuffACR(:,:,var_normc) * wMat * signal');
%             cleanDataACR = cleanDataACR(chanind(1:end-2),:)';
%             signal1ACR.signal = [cleanDataACR, force];
%             [vafACR, badLogicalACR] = cleanVAF(signal1ACR);
%             vafMatACR(var_normc,1) = mean(vafACR);
%             vafMatACR(var_normc,2) = sum(badLogicalACR);
%             vafBadLogicalACR(var_normc,:) = badLogicalACR;
%             
%             %emg
%             cleanDataEMG = signal' - (wMat_buffbuffEMG(:,:,var_normc) * wMat * signal');
%             cleanDataEMG = cleanDataEMG(chanind(1:end-2),:)';
%             signal1EMG.signal = [cleanDataEMG, force];
%             [vafEMG, badLogicalEMG] = cleanVAF(signal1EMG);
%             vafMatEMG(var_normc,1) = mean(vafEMG);
%             vafMatEMG(var_normc,2) = sum(badLogicalEMG);
%             vafBadLogicalEMG(var_normc,:) = badLogicalEMG;
        end
%         signal1.signal=signal1.signal(:,chanind);
%         signal1.parameters.ChannelNames.Value=chan;        
        
        
        
        %% SAVING
        
%         pathname = fileparts('Results/');
        
        vafname_norm=sprintf('vafNormGE_%s_%d.mat',trial(outer),ch);
        vafname_badLogical_norm=sprintf('vafNormBadLogicalGE_%s_%d.mat',trial(outer),ch);
        matfile_vafnorm=fullfile(pathname,vafname_norm);
        matfile_badlnorm=fullfile(pathname,vafname_badLogical_norm);
        save(matfile_vafnorm, 'vafMatNorm');
        save(matfile_badlnorm,'vafBadLogicalNorm' );
        
%         vafname_acr=sprintf('vafACRGE_%s_%d',trial(outer),ch);
%         vafname_badLogical_ACR=sprintf('vafACRBadLogicalGE_%s_%d.mat',trial(outer),ch);
%         matfile_vafacr=fullfile(pathname,vafname_acr);
%         matfile_badlacr=fullfile(pathname,vafname_badLogical_ACR);
%         save(matfile_vafacr, 'vafMatACR');
%         save(matfile_badlacr,'vafBadLogicalACR' );
%         
%         vafname_EMG=sprintf('vafEMGGE_%s_%d',trial(outer),ch);
%         vafname_badLogical_EMG=sprintf('vafEMGBadLogicalGE_%s_%d.mat',trial(outer),ch);
%         matfile_vafEMG=fullfile(pathname,vafname_EMG);
%         matfile_badlEMG=fullfile(pathname,vafname_badLogical_EMG);
%         save(matfile_vafEMG, 'vafMatEMG');
%         save(matfile_badlEMG,'vafBadLogicalEMG' );
       close all; 
    end
end
toc



