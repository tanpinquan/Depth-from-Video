close all

W = 10;
H = 5;
segclass = zeros(W*H,1);
pairwise = sparse(W*H,W*H);
% pairwise = zeros(50,50);

unary = ones(7,25);
[X Y] = meshgrid(1:7, 1:7);
labelcost = min(4, (X - Y).*(X - Y));

for row = 0:H-1
  for col = 0:W-1
    pixel = 1+ row*W + col;
    if row+1 < H, pairwise(pixel, 1+col+(row+1)*W) = 1; end
    if row-1 >= 0, pairwise(pixel, 1+col+(row-1)*W) = 1; end 
    if col+1 < W, pairwise(pixel, 1+(col+1)+row*W) = 1; end
    if col-1 >= 0, pairwise(pixel, 1+(col-1)+row*W) = 1; end 
    if pixel < W*H/2
      unary(:,pixel) = [0 10 10 10 10 10 10]'; 
    else
      unary(:,pixel) = [10 10 10 10 0 10 10]'; 
    end
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
