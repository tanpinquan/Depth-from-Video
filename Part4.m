close all
clc
clear
addpath('GCMex/GCMex')



nImages = 50;
param.imageWindowSize = 10;
imageStep = 1;
startImage = 90;
cameraTxt = readmatrix('images/cameras.txt');

%% Load Images
for i=1:nImages
   disp(['loading image ' num2str(startImage+(i-1)*imageStep)]);
%    disp((startImage+(i-1)*imageStep)*7+1);

   img{i} = double(imread(['images/test' num2str(startImage+(i-1)*imageStep, "%.4d") '.jpg'])); 
   imgRGB.r{i} = img{i}(:,:,1);
   imgRGB.g{i} = img{i}(:,:,2);
   imgRGB.b{i} = img{i}(:,:,3);
   
   camera{i} = cameraTxt((startImage+(i-1)*imageStep)*7+1:(startImage+(i-1)*imageStep)*7+7,1:3);
   K{i} = camera{i}(1:3,:);
   R{i} = camera{i}(4:6,:);
   T{i} = camera{i}(7,:);
   

end

[H, W, ~] = size(img{1}); 
N = H*W;

nLabels = 50;
d = 0:0.0002:0.5;
d = d(1:nLabels);

%% Compute matrices
for i = 1:nImages
    for j = 1:nImages
        mat1{i,j} = K{j}*(R{j}')*(R{i})*inv(K{i});
        mat2{i,j} = K{j}*(R{j}')*(T{i}'-T{j}');
    end
end


%% Start MRF
index_i = zeros(1,N*4 -2*W -2*H);
index_j = zeros(1,N*4 -2*W -2*H);
weights = zeros(1,N*4 -2*W -2*H);
index_count = 1;
unary = zeros(nLabels,N);


%% Cost function

param.nImages = nImages;
[param.H, param.W, ~] = size(img{1}); 

param.sigma = 1;
param.sigmaDistSq = 50;
param.nLabels = 50;
param.d = 0:0.0002:0.0002*(param.nLabels-1);
[X, Y] = meshgrid(1:nLabels, 1:nLabels);

param.lambda = 500;
param.labelcost = min(20*0.0002, abs(param.d-param.d'));

%% Perform disparity initialization

param.startImage = param.imageWindowSize+1;
param.endImage = param.nImages-param.imageWindowSize;
for i=1:param.nImages
    [disparityInit{i}] = computeDisparityInit(mat1, mat2, i, img, imgRGB, param);
end


save(['disparity_' num2str(param.nImages) '_' num2str(imageStep)], 'disparityInit')


%% Perform bundle adjustment
param.lambda = 300;
refImgIndex = 5;


for i=1:param.nImages
    disparityAdjust{i} = performBundleAdjustment(mat1, mat2, disparityInit, i, img, imgRGB, param);
end

save(['disparityAdjust_' num2str(param.nImages) '_' num2str(imageStep)], 'disparityAdjust')


