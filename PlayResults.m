clc
clear
close all


nImages = 50
windowSize = 5;
imageStep = 1;

load(['disparity_' num2str(nImages) '_' num2str(imageStep)])
load(['disparityAdjust_' num2str(nImages) '_' num2str(imageStep)])
disparityAdjustSequence = []
for i=1:nImages
    disparityInitSequence(:,:,i) = (disparityInit{i});
    disparityAdjustSequence(:,:,i) = (disparityAdjust{i});
    
end



disparityInitSequence = rescale(disparityInitSequence);
disparityAdjustSequence = rescale(disparityAdjustSequence);

implay(disparityInitSequence,10)
implay(disparityAdjustSequence,10)


