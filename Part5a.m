close all
clc
clear
addpath('GCMex/GCMex')
run('vlfeat/vlfeat-0.9.21/toolbox/vl_setup')
vl_version



nImages = 2;
param.imageWindowSize = 1;
imageStep = 20;
startImage = 0;
% cameraTxt = readmatrix('images/cameras.txt');

peakThresh = 7.65;
edgeThresh = 10;
nSubsetPoints = 5;
nFundRansacPoints = 9;
errorThresh = 1;
errorFundThresh = 0.01;

nRuns = 10000;

K = [1221.2270770	0.0000000	479.5000000	
0.0000000	1221.2270770	269.5000000	
0.0000000	0.0000000	1.0000000];
%% Load Images
for i=1:nImages
   disp(['loading image ' num2str(startImage+(i-1)*imageStep)]);
%    disp((startImage+(i-1)*imageStep)*7+1);

   img{i} = double(imread(['images/test' num2str(startImage+(i-1)*imageStep, "%.4d") '.jpg'])); 
   imgRGB.r{i} = img{i}(:,:,1);
   imgRGB.g{i} = img{i}(:,:,2);
   imgRGB.b{i} = img{i}(:,:,3);
   
%    camera{i} = cameraTxt((startImage+(i-1)*imageStep)*7+1:(startImage+(i-1)*imageStep)*7+7,1:3);
%    K{i} = camera{i}(1:3,:);
%    R{i} = camera{i}(4:6,:);
%    T{i} = camera{i}(7,:);
   
    I = single(rgb2gray(uint8(img{i}))) ;
    [sift{i}.features,sift{i}.descriptors] = vl_sift(I,'PeakThresh', peakThresh, 'EdgeThresh',edgeThresh) ;
   
    
    figure(1);
    subplot(1,nImages,i); imshow(uint8(img{i})) ;
%     h1 = vl_plotsiftdescriptor(sift{i}.descriptors(:,1:20),sift{i}.features(:,1:20)) ;
%     set(h1,'color','y') ;   
    


end

[H, W, ~] = size(img{1}); 
N = H*W;


stitchedImage = img{1};
stitchedImageRef = imref2d(size(img{1}));

points{1} = [107 306; 166 321; 444 329; 219 398; 571 373; 725 483; 782 388;151 514; 665 386];
points{2} = [216 296; 277 310; 582 326; 343 391; 675 372; 939 485; 917 391;315 508; 728 388];
figure(1); subplot(1,2,1);
hold on; plot(points{1}(:,1), points{1}(:,2),'rx')
 subplot(1,2,2);
hold on; plot(points{2}(:,1), points{2}(:,2),'rx')

direction = 1
for i=2:nImages
    [matchResult{i}.matches, matchResult{i}.scores] = matchDescriptors(sift{1}.descriptors, sift{i}.descriptors,1) ;
    
    %% Homography 
    [h{i},inlierIdx{i}] = performRansac(sift{1}.features,sift{i}.features,matchResult{i}.matches,nSubsetPoints,errorThresh,nRuns,false);
    tform = projective2d(h{i});    
    [transformedImage,transformedImageRef] = imwarp(img{i},tform);
%     sift{i}.features = applyHomographyToFeatures(h{i},sift{i}.features);
%     [stitchedImage, stitchedImageRef] = mergeImages(stitchedImage,stitchedImageRef,transformedImage,transformedImageRef,'average');
%     figure; image(stitchedImageRef.XWorldLimits, stitchedImageRef.YWorldLimits, stitchedImage) 
    
    selMatch = matchResult{i}.matches(:,inlierIdx{i});

    matchLoc{1,i} = sift{1}.features(1:2,selMatch(1,:));
    matchLoc{2,i} = sift{i}.features(1:2,selMatch(2,:));
    
     %% Fundamental
%       [F{i}, inlierIdx{i}] = performFundamentalRansac(sift{1}.features,sift{i}.features,matchResult{i}.matches,nFundRansacPoints,errorFundThresh,nRuns,false);
    F{i} = computeFundamentalMat(matchLoc{1,i}',matchLoc{2,i}');

%     F{i} = computeFundamentalMat(points{1},points{i});
    
%     [E_est{i},R_est{i}, t_est{i}] = computeEssentialMat(F{i},K{1},K{i})
%     C_est{i} = -R_est{i}'*t_est{i};
    
end
imageMerge = cat(2,uint8(img{1}),uint8(img{2}));
figure; imshow(imageMerge)

selMatch = matchResult{2}.matches(:,inlierIdx{2});
matchLoc1 = sift{1}.features(1:2,selMatch(1,:));
matchLoc2 = sift{2}.features(1:2,selMatch(2,:));
matchLoc2(1,:) = matchLoc2(1,:) + size(img{1},2);

hold on ;
l = line([matchLoc1(1,:) ; matchLoc2(1,:)], [matchLoc1(2,:) ; matchLoc2(2,:)]) ;
set(l,'linewidth', 1) ;



%% Get fundamental matrix


matchLoc1 = sift{1}.features(1:2,selMatch(1,:))';
matchLoc2 = sift{2}.features(1:2,selMatch(2,:))';

%  F{2} = estimateFundamentalMatrix(matchLoc1,matchLoc2,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',1e-4);
% 
[F{2}, inliers] = estimateFundamentalMatrix(points{1},points{2},'Method','Norm8Point')

epiLines = epipolarLine(F{2},points{1});
[E_est{2},R_est{2}, t_est{2}] = computeEssentialMat(F{2},K,K);
C_est{2} = -R_est{2}'*t_est{2};
figure;
showMatchedFeatures(uint8(img{1}),uint8(img{2}), matchLoc1, matchLoc2,'montage','PlotOptions',{'ro','go','y--'});
% figure;
% showMatchedFeatures(uint8(img{1}),uint8(img{2}), matchLoc1(inliers,:),matchLoc2(inliers,:),'montage','PlotOptions',{'ro','go','y--'});
% title('Point matches after outliers were removed');

figure;
imshow(uint8(img{2}));
title('Inliers and Epipolar Lines in Second Image'); hold on;
plot(points{2}(:,1),points{2}(:,2),'go')
epiPoints = lineToBorderPoints(epiLines,size(img{2}));

line(epiPoints(:,[1,3])',epiPoints(:,[2,4])');
truesize;
%% Test F

test1 = points{1}(1,:)';
test1(end+1) = 1
test1 = [51;463;1]
test2 = points{2}(1,:)';
test2(end+1) = 1

test2'*F{2}*test1

points{1}(:,3) = 1;
points{2}(:,3) = 1;
for i = 1:length(points{1})
    disp(points{2}(i,:)*F{2}*points{1}(i,:)')
end

l = F{2}*test1;

for i=1:nImages
    figure(100);
    subplot(1,nImages,i)
    imshow(uint8(img{i}))
end

xStart = 1;
xStep = 1;

xTest = xStart:xStep:1000*xStep;

yTest = (-l(1).*xTest - l(3))/l(2);

subplot(1,2,1);hold on;plot(test1(1),test1(2),'x')
subplot(1,2,2);hold on;plot(xTest,yTest,'x')



%% Test R and t
R_est{1} = [1 0 0; 0 1 0; 0 0 1];
t_est{1} = [0; 0; 0];
C_est{1} = -R_est{1}'*t_est{1};
for i = 1:nImages
    for j = 1:nImages
        mat1{i,j} = K*R_est{j}*(R_est{i}')*inv(K);
        mat2{i,j} = K*R_est{j}*(C_est{i}-C_est{j});
    end
end




testPoint{1} = [165; 307; 1];
% testPoint{1} = [151; 514; 1];
% testPoint{1} = [101; 170; 1];
% 
% testPoint{2} = [170; 321; 1];
% 
% testPoint{3} = [780; 386; 1];


nLabels = 100;
d = 0:0.2:nLabels*0.2;
d = d(1:nLabels);


for i=1:nImages
    figure(100);
    subplot(1,nImages,i)
    imshow(uint8(img{i}))
end
for i = 1:length(testPoint)
    for j = 1:nImages
        x{i,j} = mat1{i,j}*testPoint{i} + d.*mat2{i,j};
        x{i,j} = x{i,j}./x{i,j}(3,:);
        subplot(1,nImages,j)
        hold on;plot(x{i,j}(1,:),x{i,j}(2,:),'x')
    end
end
