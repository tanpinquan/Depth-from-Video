close all
clc
clear
addpath('GCMex/GCMex')



nImages = 2;
imageStep = 1;
startImage = 0;
cameraTxt = readmatrix('images/cameras.txt');

%% Load Images
figure(1);

img{1} = double(imread('test00.jpg'));
img{2} = double(imread('test09.jpg'));

K{1} = [1221.2270770	0.0000000	479.5000000	
0.0000000	1221.2270770	269.5000000	
0.0000000	0.0000000	1.0000000];

R{1} = [1.0000000000	0.0000000000	0.0000000000	
0.0000000000	1.0000000000	0.0000000000	
0.0000000000	0.0000000000	1.0000000000];



T{1} = [0.0000000000	0.0000000000	0.0000000000];

K{2} = [1221.2270770	0.0000000	479.5000000	
0.0000000	1221.2270770	269.5000000	
0.0000000	0.0000000	1.0000000];

R{2} = [0.9998813487	0.0148994942	0.0039106989	
-0.0148907594	0.9998865876	-0.0022532664	
-0.0039438279	0.0021947658	0.9999898146];


T{2} = [-9.9909793759	0.2451742154	0.1650832670];

for i=1:nImages
   disp(startImage+(i-1)*imageStep)
   disp((startImage+(i-1)*imageStep)*7+1)

%    img{i} = double(imread(['images/test' num2str(startImage+(i-1)*imageStep, "%.4d") '.jpg'])); 
   imgRGB.r{i} = img{i}(:,:,1);
   imgRGB.g{i} = img{i}(:,:,2);
   imgRGB.b{i} = img{i}(:,:,3);
   
%    camera{i} = cameraTxt((startImage+(i-1)*imageStep)*7+1:(startImage+(i-1)*imageStep)*7+7,1:3);
%    K{i} = camera{i}(1:3,:);
%    R{i} = camera{i}(4:6,:);
%    T{i} = camera{i}(7,:);
%    
   subplot(1,nImages,i); imshow(uint8(img{i}))
end

[H, W, ~] = size(img{1}); 
N = H*W;

%% Test points
% testPoint{1} = [165; 307; 1];
% testPoint{2} = [504; 331; 1];
% testPoint{3} = [918; 482; 1];


testPoint{1} = [641; 253; 1];
testPoint{2} = [481; 333; 1];
% testPoint{3} = [854; 389; 1];
% testPoint{4} = [806; 363; 1];

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

%% Test Point line plot
% linespec = ['rx', bx]
for i = 1:length(testPoint)
    for j = 1:nImages
        x{i,j} = mat1{i,j}*testPoint{i} + d.*mat2{i,j};
        x{i,j} = x{i,j}./x{i,j}(3,:);
        subplot(1,nImages,j)
        hold on;plot(x{i,j}(1,:),x{i,j}(2,:),'x')
    end
end



%% Start MRF
index_i = zeros(1,N*4 -2*W -2*H);
index_j = zeros(1,N*4 -2*W -2*H);
weights = zeros(1,N*4 -2*W -2*H);
index_count = 1;
unary = zeros(nLabels,N);


%% Cost function
lambda = 0.1;
sigma = 10;


param.nImages = nImages;
[param.H, param.W, ~] = size(img{1}); 

param.imageWindowSize = 1;
param.sigma = 1;
param.sigmaDistSq = 50;
param.nLabels = 50;
param.d = 0:0.0002:0.0002*(param.nLabels-1);
[X, Y] = meshgrid(1:nLabels, 1:nLabels);

param.lambda = 1250;
param.labelcost = min(20*0.0002, abs(param.d-param.d'));

for i=1:param.nImages
    [disparityInit{i}] = computeDisparityInit(mat1, mat2, i, img, imgRGB, param);
end

save(['disparity_' num2str(param.nImages) '_' num2str(imageStep)], 'disparityInit')
% load('disparityInit5.mat')

for i=1:param.nImages
% for i=param.startImage:param.endImage
    disparityAdjust{i} = performBundleAdjustment(mat1, mat2, disparityInit, i, img, imgRGB, param);
end

% refImgIndex = 1;
% disparityAdjust = performBundleAdjustment(mat1, mat2, disparityInit, refImgIndex, img, imgRGB, param);


