function [labelMat] = classifyICs(w_mix)
% labelMat is a 7 X 1 cell where every element represents the respective
% prediction
% real training data
wmat4 = load('icaw_ge404_65_115.mat'); wmat4 = wmat4.ica_weights;
wmat5 = load('icaw_ge405_65_115.mat'); wmat5 = wmat5.ica_weights;
wmat6 = load('icaw_ge406_65_115.mat'); wmat6 = wmat6.ica_weights;
wmat7 = load('icaw_ge407_65_115.mat'); wmat7 = wmat7.ica_weights;
wmat8 = load('icaw_ge408_65_115.mat'); wmat8 = wmat8.ica_weights;
wmat9 = load('icaw_ge409_65_115.mat'); wmat9 = wmat9.ica_weights;
wmat10 = load('icaw_ge410_65_115.mat'); wmat10 = wmat10.ica_weights;
% weightMat = [wmat4 wmat5 wmat6 wmat7 wmat8 wmat9 wmat10]'; % 896x128

% creating labelMat variable, 896x1620 collection of labels for each of
% 1620 classifiers across the 896 test ICs
% labelMat = zeros(896, 1800);
labelMat4 = wmat4 * w_mix';
labelMat5 = wmat5 * w_mix';
labelMat6 = wmat6 * w_mix';
labelMat7 = wmat7 * w_mix';
labelMat8 = wmat8 * w_mix';
labelMat9 = wmat9 * w_mix';
labelMat10 = wmat10 * w_mix';

labelMat4(labelMat4 < 0) = -1;
labelMat4(labelMat4 >= 0) = 1;

labelMat5(labelMat5 < 0) = -1;
labelMat5(labelMat5 >= 0) = 1;

labelMat6(labelMat6 < 0) = -1;
labelMat6(labelMat6 >= 0) = 1;

labelMat7(labelMat7 < 0) = -1;
labelMat7(labelMat7 >= 0) = 1;

labelMat8(labelMat8 < 0) = -1;
labelMat8(labelMat8 >= 0) = 1;

labelMat9(labelMat9 < 0) = -1;
labelMat9(labelMat9 >= 0) = 1;

labelMat10(labelMat10 < 0) = -1;
labelMat10(labelMat10 >= 0) = 1;


% disp(size(labelMat4))
% labelMat=cell;
labelMat{1,1}=labelMat4;
labelMat{2,1}=labelMat5;
labelMat{3,1}=labelMat6;
labelMat{4,1}=labelMat7;
labelMat{5,1}=labelMat8;
labelMat{6,1}=labelMat9;
labelMat{7,1}=labelMat10;
% labelMat=[labelMat4,labelMat5,labelMat6,labelMat7,labelMat8,labelMat9,labelMat10];


% labelMat(labelMat < 0) = -1;
% labelMat(labelMat >= 0) = 1;

end