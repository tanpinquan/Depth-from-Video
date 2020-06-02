% clc 
% clear
% close all
% 
% img1 = imread('im2.png');
% img2 = imread('im6.png');
% img1 = double(img1);
% img2 = double(img2);
% % img1 = img1(1:3:end,1:3:end,:);
% % img2 = img2(1:3:end,1:3:end,:);
% 
% img1 = img1(1:100,end-100:end,:);
% img2 = img2(1:100,end-100:end,:);
% figure; imshow(uint8([img1 img2]));
% [H,W,~] = size(img1);
% 
% 
% lambda = 10;
% dMax = floor(W/2);
% 
% nLabels = 2*W-1;
% 
% 
% N = H*W;
% segclass = zeros(N,1);
% pairwise = sparse(N,N);
% unary = zeros(nLabels,N);
% 
% [X, Y] = meshgrid(1:nLabels, 1:nLabels);
% labelcost = (X - Y).*(X - Y);
% 
% for row = 0:H-1
%     disp(row)
%     
%     for col = 0:W-1
%         pixel = 1+ row*W + col;
%         if col+1 < W, pairwise(pixel, 1+row*W+(col+1)) = lambda; end        
%         if row+1 < H, pairwise(pixel, 1+(row+1)*W+col) = lambda; end
%         if row-1 >= 0, pairwise(pixel, 1+col+(row-1)*W) = lambda; end 
%         if col-1 >= 0, pairwise(pixel, 1+(col-1)+row*W) = lambda; end 
%         
% %         a = computeDepthDataTerm(img1(row+1,col+1,:),  img1(row+1,:,:), col, W);
%         unary(:,pixel) = computeDepthDataTerm(img1(row+1,col+1,:),  img2(row+1,:,:), col, W);
%     end 
% end
% 
% [labels E Eafter] = GCMex(segclass, single(unary), pairwise, single(labelcost),0);
% 
% fprintf('E: %d (should be 260), Eafter: %d (should be 44)\n', E, Eafter);
% fprintf('unique(labels) should be [0 4] and is: [');
% fprintf('%d ', unique(labels));
% fprintf(']\n');
% 
% imageLabels = reshape(labels,[W,H]);
% imageLabels = imageLabels';
% imageLabels2 = 1./abs(imageLabels - (W+1));
% figure;imshow(uint8(imageLabels));
% 
% 



addpath('GCMex/GCMex')
close all

img = imread('bayes_in.jpg');
img = double(img);
sourceColour = [0; 0; 255]; % blue = foreground
sinkColour = [245; 210; 110]; % yellow = background
nLabels = 2;
lambda = 100;


% img = img(1:400,1:400,:);

figure; subplot(1,2,1);
imshow(uint8(img));
[H, W, ~] = size(img); 
N = H*W;

% W = 100;
% H = 100;
segclass = zeros(W*H,1);
pairwise = sparse(W*H,W*H);
% pairwise = zeros(50,50);

unary = zeros(2,25);
[X Y] = meshgrid(1:2, 1:2);
labelcost = min(4, (X - Y).*(X - Y));

for row = 0:H-1
  for col = 0:W-1
    pixel = 1+ row*W + col;
    if row+1 < H, pairwise(pixel, 1+col+(row+1)*W) = lambda; end
    if row-1 >= 0, pairwise(pixel, 1+col+(row-1)*W) = lambda; end 
    if col+1 < W, pairwise(pixel, 1+(col+1)+row*W) = lambda; end
    if col-1 >= 0, pairwise(pixel, 1+(col-1)+row*W) = lambda; end 
    
    unary(:,pixel) = computeDataTerm(img(row+1,col+1), sourceColour, sinkColour);
%     if pixel < W*H/2
%       unary(:,pixel) = [0 10]'; 
% %       unary(:,pixel) = [0 10 10 10 10 10 10]'; 
%       
%     else
%       unary(:,pixel) = [10 0]';        
% %       unary(:,pixel) = [10 10 10 10 0 10 10]'; 
%     end
  end
end

[labels E Eafter] = GCMex(segclass, single(unary), pairwise, single(labelcost),0);

fprintf('E: %d (should be 260), Eafter: %d (should be 44)\n', E, Eafter);
fprintf('unique(labels) should be [0 4] and is: [');
fprintf('%d ', unique(labels));
fprintf(']\n');

imageLabels = logical(reshape(labels,[W,H]));
imageLabels = imageLabels';
subplot(1,2,2); imshow(imageLabels);


