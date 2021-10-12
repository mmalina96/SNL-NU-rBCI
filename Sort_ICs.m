% Code to plot ICs

% SJ
% load 'ChanLocs_128.mat'
% nchan = 128;
% central = [4:6 8:11 13:15 18:21 24:26 37:40 43:44 47:50 52:54 56:59 74:79 84:89 101:106 110:115];

% load 'icaw_an_33_65_115.mat'
% load icaw_an_33_34_upto300Hz.mat
% load 'ChanLocs_128.mat'
% nchan = 128;
% central = [4:6 8:11 13:15 18:21 24:26 37:40 43:44 47:50 52:54 56:59 74:79 84:89 101:106 110:115];
% central = [4:6 8:11 13:15 18:21 24:26 42:47 52:57 69:72 75:76 79:82 84:86 88:91 101:106 110:115];

% load 'icaw_sj209_upto300Hz.mat'
% load 'ChanLocs_128.mat'
% nchan = 128;
% central = [4:6 8:11 13:15 18:21 24:26 37:40 43:44 47:50 52:54 56:59 74:79 84:89 101:106 110:115];

% load 'icaw_rs11.mat'
% load 'chanlocs_RS_updated.mat'
% nchan = 160;
% central = [42:48 60 73:77 87:88 94:103 106:114 119:120 126:133 138:141 150:156 162:167 170:175 180:184 189:197]-41;

% load 'icaw_rs11.mat'
% load 'chanlocs_RS_updated.mat'
% nchan = 160;
% central = [42:48 60 73:77 87:88 94:103 106:114 119:120 126:133 138:141 150:156 162:167 170:175 180:184 189:197]-41;

% load 'icaw_cr5_65_115.mat'
% load 'ChanLocs_128.mat'
% nchan = 128;
% central = [4:6 8:11 13:15 18:21 24:26 37:40 43:44 47:50 52:54 56:59 74:79 84:89 101:106 110:115];

% load 'icaw_fr27.mat'
% load 'ChanLocs_128.mat'
% nchan = 128;
% central = [4:6 8:11 13:15 18:21 24:26 37:40 43:44 47:50 52:54 56:59 74:79 84:89 101:106 110:115];

cd 'C:\Users\micha\Documents\SNL_NU\FromMukta'
wMat = load('icaw_ge404_65_115.mat');
signal1 = load('GE_handS004R04_mvcut.mat');
force = signal1.signal(:,129:130);
ica_weights = wMat.ica_weights;
% % load icaw_cm_703_under300.mat
load 'ChanLocs_128.mat'
nchan = 128;
central = [4:6 8:11 13:15 18:21 24:26 37:40 43:44 47:50 52:54 56:59 74:79 84:89 101:106 110:115];
% 
% load icaw_ge410_fRobert.mat
% nchan = 108;
% % central = [10:36 38:42 44:55 57:62 73 77 79:83]; %GE404_freqRob
% central = [11:33 36 38:42 44:54 56:60 70 76 78:82 99]; %GE410_freqRob

% load 'icaw_mb218_65_115.mat'
% load 'eeglab_locs_mb2.mat'
% locs = eeglab_locs;
% nchan=64;
% central = [4:6 8:11 13:15 18:21 24:26 37:40 43:44 47:50 52:54 56:59 62];


%% ACR
nICs = size(ica_weights,1);
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

periph = setdiff(1:nchan,central);
keep = [];
throw = [];
for ic =1:size(ica_weights,1)
    mean_central = mean(abs(ica_weights(ic, central)));
    mean_periph = mean(abs(ica_weights(ic,periph)));
    if mean_central>(mean_periph*1)
        keep = [keep ic];
    else
        throw = [throw ic];
    end
end

resultMat = zeros(128,1);
resultMat(keep) = 1;
resultMat(throw) = -1;




hemi={'Fz','F4','F6','FFC4h','FFC6h','FC2','FC4','FC6','FCC4h','FCC6h','C4','C6','CCP4h','CCP6h','CP4','CP6','force','force2'};
homoLog = {'F1', 'F3', 'F5', 'FFC3h', 'FFC5h', 'FC1', 'FC3', 'FC5', 'FCC3h', 'FCC5h', 'C3', 'C5', 'CCP4h', 'CCP5h', 'CP3', 'CP5', 'force', 'force2'}; 

chan = hemi;

chancount=1;
tmp=signal1.parameters.ChannelNames.Value;
for count =1:130
    namecell=tmp(count,1);
    st=strcmp(namecell{1,1},chan); %change depending on homoLog or homoLog
    if sum(st)==1
        chanind(chancount)=count;
        chancount=chancount+1;
    end
end


dirtyDataACRsignal = signal1.signal(:,1:128);
for i = 1:128
    if resultMat(i) == -1
        wMatClean(i,:) = zeros(128,1);
    else
        wMatClean(i,:) = wMat.ica_weights(i,:);
    end
end
cleanDataACR = dirtyDataACRsignal' - (wMatClean * wMat.ica_weights * dirtyDataACRsignal');
cleanDataACR = cleanDataACR(chanind(1:end-2),:)'; %chanind from cleanICs
signal1.signal = cleanDataACR;
signal1.parameters.ChannelNames.Value = chan;


signal1.signal = [signal1.signal force];

vaf = cleanVAF(signal1)



% %% plotting
% for ic = keep(1:5)
%     figure
%     scatter3(allx,ally,allz,100,ica_weights(ic,:),'filled')
%     for loc = 1:128
% %         RS = [[1 19 23 64 71 100 104 132 137]]; CR = [5 14 25 30 62 65 13 15];
% %         MB2 = [2 6 16 26 29 31 49 15 17];
%         lab = locs(loc).labels;
%         text(allx(1,loc),ally(1,loc),allz(1,loc),lab,'HorizontalAlignment','center','FontWeight','bold')
%     hold on
%     end
%     colormap(jet)
%     axis off
%     maxval = max(abs(min(min(ica_weights(ic,:)))), abs(max(max(ica_weights(ic,:)))));
%     caxis([maxval*-1 maxval])
%     colorbar()
%     title(['Keep: IC#' num2str(ic)])
% end
% 
% for ic = throw(1:5)
%     figure
%     scatter3(allx,ally,allz,100,ica_weights(ic,:),'filled')
%     for loc = 1:128
% %         RS = [[1 19 23 64 71 100 104 132 137]]; CR = [5 14 25 30 62 65 13 15];
% %         MB2 = [2 6 16 26 29 31 49 15 17];
%         lab = locs(loc).labels;
%         text(allx(1,loc),ally(1,loc),allz(1,loc),lab,'HorizontalAlignment','center','FontWeight','bold')
%     hold on
%     end
%     colormap(jet)
%     axis off
%     maxval = max(abs(min(min(ica_weights(ic,:)))), abs(max(max(ica_weights(ic,:)))));
%     caxis([maxval*-1 maxval])
%     colorbar()
%     title(['Throw: IC#' num2str(ic)])
% end

    