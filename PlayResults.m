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

%     if(i>windowSize && i <nImages-windowSize+1)
%         disparityAdjustSequence(:,:,end+1) = (disparityAdjust{i});
%     end
    
end



disparityInitSequence = rescale(disparityInitSequence);
disparityAdjustSequence = rescale(disparityAdjustSequence);
implay(disparityInitSequence)
implay(disparityAdjustSequence)

% outputVideo = VideoWriter(fullfile('disparity_adjust2.avi'));
% outputVideo.FrameRate = 10;
% open(outputVideo)
% 
% for i = 10:30
%    writeVideo(outputVideo,disparityAdjustSequence(:,:,i))
% end
% 
% close(outputVideo)

