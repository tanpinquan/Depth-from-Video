close all
clc
clear
addpath('GCMex/GCMex')

img1 = imread('test00.jpg');
img2 = imread('test09.jpg');
img1 = double(img1);
img2 = double(img2);


img2RGB.r = img2(:,:,1);
img2RGB.g = img2(:,:,2);
img2RGB.b = img2(:,:,3);


[H, W, ~] = size(img1); 
N = H*W;

%% K, R and T
K = [1221.2270770	0.0000000	479.5000000	
0.0000000	1221.2270770	269.5000000	
0.0000000	0.0000000	1.0000000];

R = [1.0000000000	0.0000000000	0.0000000000	
0.0000000000	1.0000000000	0.0000000000	
0.0000000000	0.0000000000	1.0000000000];

T = [0.0000000000	0.0000000000	0.0000000000];

Kprime = [1221.2270770	0.0000000	479.5000000	
0.0000000	1221.2270770	269.5000000	
0.0000000	0.0000000	1.0000000];

Rprime = [0.9998813487	0.0148994942	0.0039106989	
-0.0148907594	0.9998865876	-0.0022532664	
-0.0039438279	0.0021947658	0.9999898146];

Tprime = [-9.9909793759	0.2451742154	0.1650832670];

%% Swap K, R and T
% img1 = double(img1);
imgTemp = img2;

img2 = img1;
img1 = imgTemp;

img2RGB.r = img2(:,:,1);
img2RGB.g = img2(:,:,2);
img2RGB.b = img2(:,:,3);

Kprime = [1221.2270770	0.0000000	479.5000000	
0.0000000	1221.2270770	269.5000000	
0.0000000	0.0000000	1.0000000];

Rprime = [1.0000000000	0.0000000000	0.0000000000	
0.0000000000	1.0000000000	0.0000000000	
0.0000000000	0.0000000000	1.0000000000];

Tprime = [0.0000000000	0.0000000000	0.0000000000];

K = [1221.2270770	0.0000000	479.5000000	
0.0000000	1221.2270770	269.5000000	
0.0000000	0.0000000	1.0000000];

R = [0.9998813487	0.0148994942	0.0039106989	
-0.0148907594	0.9998865876	-0.0022532664	
-0.0039438279	0.0021947658	0.9999898146];

T = [-9.9909793759	0.2451742154	0.1650832670];

nLabels = 100;

[X Y] = meshgrid(1:nLabels, 1:nLabels);

%% Test Point
testPoint = [725; 483; 1];

d = 0:0.0001:0.5;
d = d(1:nLabels);
mat1 = Kprime*(Rprime')*(R)*inv(K)
mat2 = Kprime*(Rprime')*(T'-Tprime')

% lambda2= 0.2;
% labelcost2 = min(9, (X - Y).*(X - Y));
lambda = 2e7;
labelcost = min(9*0.0001^2, (d-d').*(d-d'));

for i = 1:length(d)
    testPointPrime(:,i) = mat1*testPoint+d(i)*mat2;
end


testPointPrime = testPointPrime./testPointPrime(3,:);

validPoints = testPointPrime(1,:)>1 & testPointPrime(1,:)<W & testPointPrime(2,:)>1 & testPointPrime(2,:)<H;

validTestPointPrime = testPointPrime(:,validPoints);
validTestPointPrimeRounded = round(validTestPointPrime);
% testPointPrime(validPoints,:);
figure; 
subplot(1,2,1);
imshow(uint8([img1]));
subplot(1,2,2);
imshow(uint8(img2));

subplot(1,2,1)
hold on;plot(testPoint(1),testPoint(2),'rx', 'markersize', 8);
subplot(1,2,2)
hold on;plot(validTestPointPrime(1,:),validTestPointPrime(2,:),'r.');

testDataTerm = computeUnrectifiedDepthDataTerm(testPoint, d, mat1, mat2, img1, img2RGB, W, H, nLabels, 'norm', 10);

priorCost = testDataTerm+lambda*labelcost(:,1);
ind = round(mean(find(priorCost == min(priorCost))));

validTestPointPrime(:,ind)

%% Start MRF
index_i = zeros(1,N*4 -2*W -2*H);
index_j = zeros(1,N*4 -2*W -2*H);
weights = zeros(1,N*4 -2*W -2*H);
index_count = 1;


unary = zeros(nLabels,N);

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
    
    x = [col+1; row+1;1];
    unary(:,pixel) = computeUnrectifiedDepthDataTerm(x, d, mat1, mat2, img1, img2RGB, W, H, nLabels, 'norm', 10);
%     unary(:,pixel) = computeDepthDataTerm(img1(row+1,col+1,:),  img2(row+1,:,:), col, nLabels);
    
%     priorCost = unary(:,pixel)+lambda*labelcost(:,1);
%     index = round(mean(find(priorCost == min(priorCost))));
    
%     [val index] = min(unary(:,pixel)+lambda*labelcost(:,1));
    
    [val index] = min(unary(:,pixel));    
    segclass(pixel) = index-1;

  end
end

priorProb = reshape(segclass,[W,H]);
priorProb = priorProb';
figure(2);subplot(1,2,1); imshow(rescale(priorProb));
title('Data term')


pairwise2 = sparse(index_i, index_j, weights);

dataTerm = unary';

disp(['graphcuts start' datestr(now, 'dd/mm/yy-HH:MM')])

[labels E Eafter] = GCMex(segclass, single(unary), pairwise2, single(labelcost),0);

disp(['graphcuts end' datestr(now, 'dd/mm/yy-HH:MM')])

fprintf('E: %d (should be 260), Eafter: %d (should be 44)\n', E, Eafter);
fprintf('unique(labels) should be [0 4] and is: [');
fprintf('%d ', unique(labels));
fprintf(']\n')

imageLabels = (reshape(labels,[W,H]));
imageLabels = imageLabels';
figure(2); subplot(1,2,2); imshow(rescale(imageLabels));
title('Depth term')


