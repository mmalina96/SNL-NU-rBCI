% we wanto define characteristics of IC weight matrices, independent of
% location of activation or whether it is EMG/EEG. Just want to find
% typical means, SDs, maxes and mins of dataset for future modeling. 
% Have 896 ICs to sample from (7 x 128)
%function [acrMixMat, emgMixMat] =  testDataLoad 
cd C:\Users\micha\Documents\SNL_NU\FromMukta
load 'ChanLocs_128.mat'

cd C:\Users\micha\Documents\SNL_NU\FromMukta\intentEMGs

%loading in all (that I currently have) datasets
wmat4 = load('ica_intent_CM_21_65_115.mat'); wmat4 = wmat4.ica_weights;
wmat5 = load('ica_intent_CM_31_65_115.mat'); wmat5 = wmat5.ica_weights;
wmat6 = load('ica_intent_CM_42_65_115.mat'); wmat6 = wmat6.ica_weights;
wmat7 = load('ica_intent_CM_51_65_115.mat'); wmat7 = wmat7.ica_weights;
wmat8 = load('ica_intent_CM_61_65_115.mat'); wmat8 = wmat8.ica_weights;
wmat9 = load('ica_intent_CS_65_115.mat'); wmat9 = wmat9.ica_weights;

%putting them all into one 3D mat
weights(:,:,1) = wmat4;
weights(:,:,2) = wmat5;
weights(:,:,3) = wmat6;
weights(:,:,4) = wmat7;
weights(:,:,5) = wmat8;
weights(:,:,6) = wmat9;

weightMat = [wmat4 wmat5 wmat6 wmat7 wmat8 wmat9]';
nchan = 128;

%% ACR
nICs = 128;
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

keep = [];
throw = [];
for ic =1:size(weightMat,1)
    mean_central = mean(abs(weightMat(ic, central)));
    mean_periph = mean(abs(weightMat(ic,periph)));
    if mean_central>(mean_periph*1)
        keep = [keep ic];
    else
        throw = [throw ic];
    end
end

acrMixMat = zeros(768,129);
acrMixMat(keep,1:128) = weightMat(keep,:);
acrMixMat(keep,129) = ones(length(keep),1);
acrMixMat(throw,1:128) = weightMat(throw,:);
acrMixMat(throw,129) = -1 .* ones(length(throw),1);

emgMixMat = zeros(768,129);
emgMixMat(:,1:128) = weightMat;
emgMixMat(:,129) = -1 .* ones(length(weightMat),1);

save('acrMixMat.mat', 'acrMixMat')
save('emgMixMat.mat', 'emgMixMat')

% keepL = ones(length(keep),1);
% throwL = -1 .* ones(length(throw),1);
% 
% testMat = [keep' keepL;  throw' throwL];
% 
% testData = zeros(896, 129);
% testData(1:size(weightMat,1),1:128) = weightMat;
% for i = 1:size(weightMat,1)
%     testData(testMat(i,1),129) = testMat(i,2);
% end
