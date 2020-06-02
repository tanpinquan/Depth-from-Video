close all
clc
clear
addpath('GCMex/GCMex')

img1 = imread('im6.png');
img2 = imread('im2.png');
img1 = double(img1);
img2 = double(img2);

figure; 
subplot(1,2,1);
imshow(uint8([img1]));
subplot(1,2,2);
imshow(uint8(img2));
[H, W, ~] = size(img1); 
N = H*W;


segclass = zeros(W*H,1);

index_i = zeros(1,N*4 -2*W -2*H);
index_j = zeros(1,N*4 -2*W -2*H);
weights = zeros(1,N*4 -2*W -2*H);
index_count = 1;

maxDisplacement = 60;
nLabels = 2*maxDisplacement+1;

[X Y] = meshgrid(1:nLabels, 1:nLabels);
maxCost = 6.5
labelcost = min(maxCost, (X - Y).*(X - Y));


method = 'abs'
lambda = 5;
sigma = 10;



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
    
  
    unary(:,pixel) = computeDepthDataTerm(img1(row+1,col+1,:),  img2(row+1,:,:), col, maxDisplacement, method, sigma);

    [val index] = min(unary(:,pixel));
    segclass(pixel) = index-1;

  end
end

priorProb = (segclass - (maxDisplacement+1));
priorProb = abs(priorProb);
priorProb = reshape(priorProb,[W,H]);
priorProb = priorProb';

pairwise2 = sparse(index_i, index_j, weights);

dataTerm = unary';

disp(['Starting graphcuts: ' datestr(now, 'dd/mm/yy-HH:MM')])

[labels E Eafter] = GCMex(segclass, single(unary), pairwise2, single(labelcost),0);

fprintf('E: %d , Eafter: %d \n', E, Eafter);
fprintf('%d ', unique(labels));
fprintf(']\n');

imageLabels = (reshape(labels,[W,H]));
imageLabels = imageLabels - (maxDisplacement+1);
imageLabels = abs(imageLabels);
imageLabels = imageLabels';

figure;
imshow(rescale(priorProb)); title('Data Term')

figure;
imshow(rescale(imageLabels)); title('Disparity map')

disp(['Complete: ' datestr(now, 'dd/mm/yy-HH:MM')])
