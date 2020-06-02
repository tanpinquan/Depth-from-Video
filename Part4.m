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
   disp(startImage+(i-1)*imageStep)
   disp((startImage+(i-1)*imageStep)*7+1)

   img{i} = double(imread(['images/test' num2str(startImage+(i-1)*imageStep, "%.4d") '.jpg'])); 
   imgRGB.r{i} = img{i}(:,:,1);
   imgRGB.g{i} = img{i}(:,:,2);
   imgRGB.b{i} = img{i}(:,:,3);
   
   camera{i} = cameraTxt((startImage+(i-1)*imageStep)*7+1:(startImage+(i-1)*imageStep)*7+7,1:3);
   K{i} = camera{i}(1:3,:);
   R{i} = camera{i}(4:6,:);
   T{i} = camera{i}(7,:);
   
   % figure(1);
%    subplot(1,nImages,i); imshow(uint8(img{i}))
end

[H, W, ~] = size(img{1}); 
N = H*W;

%% Test points
testPoint{1} = [165; 307; 1];
testPoint{2} = [504; 331; 1];
testPoint{3} = [918; 482; 1];


testPoint{1} = [641; 253; 1];
testPoint{2} = [481; 333; 1];
testPoint{3} = [854; 389; 1];
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
% for i = 1:length(testPoint)
%     for j = 1:nImages
%         x{i,j} = mat1{i,j}*testPoint{i} + d.*mat2{i,j};
%         x{i,j} = x{i,j}./x{i,j}(3,:);
%         subplot(1,nImages,j)
%         hold on;plot(x{i,j}(1,:),x{i,j}(2,:),'x')
%     end
% end



%% Start MRF
index_i = zeros(1,N*4 -2*W -2*H);
index_j = zeros(1,N*4 -2*W -2*H);
weights = zeros(1,N*4 -2*W -2*H);
index_count = 1;
unary = zeros(nLabels,N);


%% Cost function
% lambda = 0.1;
% sigma = 10;

% labelcost = min(9, (X - Y).*(X - Y));

param.nImages = nImages;
[param.H, param.W, ~] = size(img{1}); 

% param.lambda = 500;
param.sigma = 1;
param.sigmaDistSq = 50;
param.nLabels = 50;
param.d = 0:0.0002:0.0002*(param.nLabels-1);
[X, Y] = meshgrid(1:nLabels, 1:nLabels);
% param.lambda = 0.1;
% param.labelcost = min(9, (X - Y).*(X - Y));
param.lambda = 500;
param.labelcost = min(20*0.0002, abs(param.d-param.d'));

%% Perform disparity initialization

param.startImage = param.imageWindowSize+1;
param.endImage = param.nImages-param.imageWindowSize;
for i=1:param.nImages
    [disparityInit{i}] = computeDisparityInit(mat1, mat2, i, img, imgRGB, param);
end


save(['disparity_' num2str(param.nImages) '_' num2str(imageStep)], 'disparityInit')
% load(['disparity_' num2str(nImages) '_' num2str(imageStep)])


%% Perform bundle adjustment
param.lambda = 300;
refImgIndex = 5;

% disparityAdjust = performBundleAdjustment(mat1, mat2, disparityInit, refImgIndex, img, imgRGB, param);

for i=1:param.nImages
% for i=param.startImage:param.endImage
    disparityAdjust{i} = performBundleAdjustment(mat1, mat2, disparityInit, i, img, imgRGB, param);
end

save(['disparityAdjust_' num2str(param.nImages) '_' num2str(imageStep)], 'disparityAdjust')

% depthMap = 1./disparityAdjust
% figure(310+refImgIndex+1);surf(flipud(-depthMap),'EdgeColor','None');
% colormap('gray');








% 
% for row = 0:H-1
%     disp(row)
%     for col = 0:W-1  
% %     disp(col)
% 
%     pixel = 1+ row*W + col;
%     
%     
%     if row+1 < H 
%         index_i(index_count) = pixel;
%         index_j(index_count) = 1+col+(row+1)*W;
%         weights(index_count) = lambda;
%         index_count = index_count+1;
%     end
%     if row-1 >= 0 
%         index_i(index_count) = pixel;
%         index_j(index_count) = 1+col+(row-1)*W;
%         weights(index_count) = lambda;
%         index_count = index_count+1;
%         
%     end 
%     if col+1 < W 
%         index_i(index_count) = pixel;
%         index_j(index_count) = 1+(col+1)+row*W;
%         weights(index_count) = lambda;   
%         index_count = index_count+1;
% 
%     end
%     if col-1 >= 0
%         index_i(index_count) = pixel;
%         index_j(index_count) = 1+(col-1)+row*W;
%         weights(index_count) = lambda; 
%         index_count = index_count+1;
%         
%     end    
%     
%     x = [col+1; row+1;1];
% %     unary(:,pixel) = computeVideoDepthDataTerm(x, d, mat1, mat2, refImgIndex, img, imgRGB, W, H, nImages, nLabels, sigma);
% %     [val index] = min(unary(:,pixel));    
% %     segclass(pixel) = index-1;
%     
%     unary(:,pixel) = computeVideoDepthDataTerm2(x, mat1, mat2, refImgIndex, img, imgRGB,param);
%     [val index] = min(unary(:,pixel));    
%     segclass(pixel) = index-1;
%     end  
% end
% priorProb = reshape(segclass,[W,H]);
% priorProb = priorProb';
% figure(2);subplot(1,2,1); imshow(rescale(priorProb));
% 
% pairwise2 = sparse(index_i, index_j, weights);
% 
% dataTerm = unary';
% 
% disp(datestr(now, 'dd/mm/yy-HH:MM'))
% 
% [labels E Eafter] = GCMex(segclass, single(unary), pairwise2, single(param.labelcost),0);
% 
% disp(datestr(now, 'dd/mm/yy-HH:MM'))
% 
% fprintf('E: %d (should be 260), Eafter: %d (should be 44)\n', E, Eafter);
% fprintf('unique(labels) should be [0 4] and is: [');
% fprintf('%d ', unique(labels));
% fprintf(']\n')
% 
% imageLabels = (reshape(labels,[W,H]));
% imageLabels = imageLabels';
% figure(2); subplot(1,2,2); imshow(rescale(imageLabels));
