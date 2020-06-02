close all
clc
clear
addpath('GCMex/GCMex')

img = imread('bayes_in.jpg');
img = double(img);
sourceColour = [0; 0; 255]; % blue = foreground
sinkColour = [245; 210; 110]; % yellow = background
nLabels = 2;
lambda = 1000;

% img = img(1:2:end, 1:2:end,:);
% img = img(1:400,1:401,:);

% imshow(uint8(img));
[H, W, ~] = size(img); 
N = H*W;

% W = 100;
% H = 100 ;
segclass = zeros(W*H,1);
% segclass = round(rand(W*H,1));

pairwise = sparse(W*H,W*H);
% pairwise = zeros(50,50);
index_i = zeros(1,N*4 -2*W -2*H);
index_j = zeros(1,N*4 -2*W -2*H);
weights = zeros(1,N*4 -2*W -2*H);
index_count = 1;

% unary = zeros(2,25);
[X Y] = meshgrid(1:2, 1:2);
labelcost = (X - Y).*(X - Y);


for row = 0:H-1
    disp(['Processing row ' num2str(row)])
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
    
    unary(:,pixel) = computeDataTerm(img(row+1,col+1,:), sourceColour, sinkColour);
    [val index] = min(unary(:,pixel));
    segclass(pixel) = index-1;
    
%     unary2(:,pixel) = computeDataTerm(img(row+1,col+1), sourceColour, sinkColour);
    
%     if(row>399)
%         unary(:,pixel) = [10 0];
% %         replacementIndex = pixel - ( 400*W ) 
% %         unary(:,pixel) = unary(:,replacementIndex);        
%     end

%     if pixel < W*H/2
%       unary(:,pixel) = [0 10 ]'; 
%     else
%       unary(:,pixel) = [10 0]'; 
%     end
  end
end

pairwise2 = sparse(index_i, index_j, weights);

dataTerm = unary';
[labels E Eafter] = GCMex(segclass, single(unary), pairwise2, single(labelcost),1);

fprintf('E: %d (should be 260), Eafter: %d (should be 44)\n', E, Eafter);
fprintf('unique(labels) should be [0 4] and is: [');
fprintf('%d ', unique(labels));
fprintf(']\n');

imageLabels = logical(reshape(labels,[W,H]));
imageLabels = imageLabels';
% 
% figure;
% imshow(imageLabels);
% set(gcf, 'Position',  [200, 200, 400, 300])
% title(['lambda = ' num2str(lambda)])

imgOut.r = 0*uint8(ones(H,W));
imgOut.g = 0*uint8(ones(H,W));
imgOut.b = 255*uint8(ones(H,W));
imgOut.r(imageLabels) = 245;
imgOut.g(imageLabels) = 210;
imgOut.b(imageLabels) = 110;

imgOut.rgb(:,:,1) = imgOut.r;
imgOut.rgb(:,:,2) = imgOut.g;
imgOut.rgb(:,:,3) = imgOut.b;
figure;imshow(imgOut.rgb);
set(gcf, 'Position',  [200, 200, 400, 300])
title(['lambda = ' num2str(lambda)])
