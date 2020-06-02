% clc 
% clear
% close all
% addpath('GCMex/GCMex')
% 
% img = imread('bayes_in.jpg');
% img = double(img);
% sourceColour = [0; 0; 255]; % blue = foreground
% sinkColour = [245; 210; 110]; % yellow = background
% nLabels = 2;
% lambda = 100;
% 
% % img = img(1:3:end,1:3:end,:);
% img = img(1:400,1:400,:);
% 
% figure; subplot(1,2,1);
% imshow(uint8(img));
% [H, W, ~] = size(img); 
% N = H*W;
% 
% % segclass = zeros(N,1);
% % pairwise = sparse(N,N);
% % unary = zeros(nLabels,N);
% % labelcost = [0 1; 
% %             1 0];
% 
% segclass = zeros(W*H,1);
% pairwise = sparse(W*H,W*H);        
% unary = zeros(2,25);
% [X Y] = meshgrid(1:2, 1:2);
% labelcost = min(4, (X - Y).*(X - Y));
% 
% for row = 0:H-1
%     for col = 0:W-1
%         pixel = 1+ row*W + col;
%         if row+1 < H, pairwise(pixel, 1+(row+1)*W+col) = lambda; end
%         if row-1 >= 0, pairwise(pixel, 1+col+(row-1)*W) = lambda; end 
%         if col+1 < W, pairwise(pixel, 1+row*W+(col+1)) = lambda; end                
%         if col-1 >= 0, pairwise(pixel, 1+(col-1)+row*W) = lambda; end 
%         
% %         a  =computeDataTerm(img(row+1,col+1,:),  sourceColour, sinkColour);
%         unary(:,pixel) = computeDataTerm(img(row+1,col+1,:),  sourceColour, sinkColour);
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
% imageLabels = logical(reshape(labels,[W,H]));
% imageLabels = imageLabels';
% subplot(1,2,2); imshow(imageLabels);
% 
% 
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


img = img(1:400,1:400,:);

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
    disp(row)
  for col = 0:W-1
    pixel = 1+ row*W + col;
    if row+1 < H, pairwise(pixel, 1+col+(row+1)*W) = lambda; end
    if row-1 >= 0, pairwise(pixel, 1+col+(row-1)*W) = lambda; end 
    if col+1 < W, pairwise(pixel, 1+(col+1)+row*W) = lambda; end
    if col-1 >= 0, pairwise(pixel, 1+(col-1)+row*W) = lambda; end 
    
    unary(:,pixel) = computeDataTerm(img(row+1,col+1), sourceColour, sinkColour);

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

