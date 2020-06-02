% close all
clc
clear
addpath('GCMex/GCMex')

img1 = imread('im2.png');
img2 = imread('im6.png');
img1 = double(img1);
img2 = double(img2);

% img1 = img1(1:2:end, 1:2:end,:);
% img2 = img2(1:2:end, 1:2:end,:);

% img = img(1:400,1:401,:);

figure; 
subplot(1,3,1);
imshow(uint8([img1]));
subplot(1,3,2);
imshow(uint8(img2));
[H, W, ~] = size(img1); 
N = H*W;

% W = 100;
% H = 100 ;
segclass = zeros(W*H,1);
% pairwise = sparse(W*H,W*H);
% pairwise = zeros(50,50);
index_i = zeros(1,N*4 -2*W -2*H);
index_j = zeros(1,N*4 -2*W -2*H);
weights = zeros(1,N*4 -2*W -2*H);
index_count = 1;

% unary = zeros(2,25);
% nLabels = 2;
maxDisplacement = 60;
nLabels = 2*maxDisplacement+1;

[X Y] = meshgrid(1:nLabels, 1:nLabels);
maxCost = 6.5
labelcost = min(maxCost, (X - Y).*(X - Y));
% labelcost = min(10, abs(X - Y));

% method = 'norm';
% lambda = 0.2;
% sigma = 10;

method = 'abs'
lambda = 5;
sigma = 10;

% lambda = 0.5;
% labelcost = abs(X - Y);


for row = 0:H-1
    disp(row)
  for col = 0:W-1
    pixel = 1+ row*W + col;

    
    if row+1 < H 
        index_i(index_count) = pixel;
        index_j(index_count) = 1+col+(row+1)*W;
        weights(index_count) = lambda;
        index_count = index_count+1;
    end
    if row-1 >= 0 
        index_i(index_count) = pixel;
        index_j(index_count) = 1+col+(row-1)*W;
        weights(index_count) = lambda;
        index_count = index_count+1;
        
    end 
    if col+1 < W 
        index_i(index_count) = pixel;
        index_j(index_count) = 1+(col+1)+row*W;
        weights(index_count) = lambda;   
        index_count = index_count+1;

    end
    if col-1 >= 0
        index_i(index_count) = pixel;
        index_j(index_count) = 1+(col-1)+row*W;
        weights(index_count) = lambda; 
        index_count = index_count+1;
        
    end
    
    
%     if row+1 < H, pairwise(pixel, 1+col+(row+1)*W) = lambda; end
%     if row-1 >= 0, pairwise(pixel, 1+col+(row-1)*W) = lambda; end 
%     if col+1 < W, pairwise(pixel, 1+(col+1)+row*W) = lambda; end
%     if col-1 >= 0, pairwise(pixel, 1+(col-1)+row*W) = lambda; end 
        
    unary(:,pixel) = computeDepthDataTerm(img1(row+1,col+1,:),  img2(row+1,:,:), col, maxDisplacement, method, sigma);

    [val index] = min(unary(:,pixel));
    segclass(pixel) = index-1;

  end
end

% priorProb = abs(segclass - (W+1));
priorProb = (segclass - (maxDisplacement+1));
priorProb = abs(priorProb);
priorProb = reshape(priorProb,[W,H]);
priorProb = priorProb';

pairwise2 = sparse(index_i, index_j, weights);

dataTerm = unary';

disp(datestr(now, 'dd/mm/yy-HH:MM'))

[labels E Eafter] = GCMex(segclass, single(unary), pairwise2, single(labelcost),0);

fprintf('E: %d (should be 260), Eafter: %d (should be 44)\n', E, Eafter);
fprintf('unique(labels) should be [0 4] and is: [');
fprintf('%d ', unique(labels));
fprintf(']\n');

imageLabels = (reshape(labels,[W,H]));
imageLabels = imageLabels - (maxDisplacement+1);
imageLabels = abs(imageLabels);
imageLabels = imageLabels';

figure;
% subplot(1,2,1); imshow(rescale(priorProb));
% subplot(1,2,2); 
imshow(rescale(imageLabels));
% title([method ', lambda: ' num2str(lambda) ' maxcost: ' num2str(maxCost)])

% figure; imagesc(1./imageLabels);

% figure;surf(flipud(imageLabels),'EdgeColor','None');
% view(2);
disp(datestr(now, 'dd/mm/yy-HH:MM'))
