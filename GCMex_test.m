close all
clear
clc

W = 5;
H = 10;
N = W*H;

lambda = 1;

segclass = zeros(W*H,1);
pairwise = sparse(W*H,W*H);
% pairwise = zeros(50,50);
index_i = zeros(1,N*4 -2*W -2*H);
index_j = zeros(1,N*4 -2*W -2*H);
weights = zeros(1,N*4 -2*W -2*H);
index_count = 1;

unary = ones(2,25);
[X Y] = meshgrid(1:2, 1:2);
labelcost = min(4, (X - Y).*(X - Y));

for row = 0:H-1
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
    
%     if row+1 < H, pairwise(pixel, 1+col+(row+1)*W) = 1; end
%     if row-1 >= 0, pairwise(pixel, 1+col+(row-1)*W) = 1; end 
%     if col+1 < W, pairwise(pixel, 1+(col+1)+row*W) = 1; end
%     if col-1 >= 0, pairwise(pixel, 1+(col-1)+row*W) = 1; end 


    
    if pixel < W*H/2
      unary(:,pixel) = [0 10 ]'; 
    else
      unary(:,pixel) = [10 0]'; 
    end
  end
end
pairwise2 = sparse(index_i, index_j, weights);

[labels E Eafter] = GCMex(segclass, single(unary), pairwise2, single(labelcost),0);

fprintf('E: %d (should be 260), Eafter: %d (should be 44)\n', E, Eafter);
fprintf('unique(labels) should be [0 4] and is: [');
fprintf('%d ', unique(labels));
fprintf(']\n');

imageLabels = logical(reshape(labels,[W,H]));
imageLabels = imageLabels';
subplot(1,2,2); imshow(imageLabels);
