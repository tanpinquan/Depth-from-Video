function [disparityLabels] = computeDisparityInit(mat1, mat2, refImgIndex, img, imgRGB, param)

%% Start MRF
N = param.W*param.H;
index_i = zeros(1,N*4 -2*param.W -2*param.H);
index_j = zeros(1,N*4 -2*param.W -2*param.H);
weights = zeros(1,N*4 -2*param.W -2*param.H);
index_count = 1;
unary = zeros(param.nLabels,N);
segclass = zeros(1,N);
%% Cost function

disp(['Data Term ' num2str(refImgIndex) ' start:' datestr(now, 'dd/mm/yy-HH:MM')])

for row = 0:param.H-1
% for row = 259:param.H-1
% for row = 329:param.H-1
    
%     disp(['Disparity init image:' num2str(refImgIndex) ', row:' num2str(row)])
    for col = 0:param.W-1  
%     for col = 699:param.W-1  
%     for col = 613:param.W-1  
        
%     disp(col)

    pixel = 1+ row*param.W + col;
    
    
    if row+1 < param.H 
        index_i(index_count) = pixel;
        index_j(index_count) = 1+col+(row+1)*param.W;
        weights(index_count) = param.lambda;
        index_count = index_count+1;
    end
    if row-1 >= 0 
        index_i(index_count) = pixel;
        index_j(index_count) = 1+col+(row-1)*param.W;
        weights(index_count) = param.lambda;
        index_count = index_count+1;
        
    end 
    if col+1 < param.W 
        index_i(index_count) = pixel;
        index_j(index_count) = 1+(col+1)+row*param.W;
        weights(index_count) = param.lambda;   
        index_count = index_count+1;

    end
    if col-1 >= 0
        index_i(index_count) = pixel;
        index_j(index_count) = 1+(col-1)+row*param.W;
        weights(index_count) = param.lambda; 
        index_count = index_count+1;
        
    end    
    
    x = [col+1; row+1;1];
    unary(:,pixel) = computeVideoDepthDataTerm(x, mat1, mat2, refImgIndex, img, imgRGB,param);
    [val index] = min(unary(:,pixel));    
    segclass(pixel) = index-1;
    
    end  
end
priorProb = reshape(segclass,[param.W,param.H]);
priorProb = priorProb';
figure(100 + refImgIndex);subplot(1,2,1); imshow(rescale(priorProb)); title('Data Term')


pairwise2 = sparse(index_i, index_j, weights);


disp(['Graph cut start:' datestr(now, 'dd/mm/yy-HH:MM')])

[labels, E, Eafter] = GCMex(segclass, single(unary), pairwise2, single(param.labelcost),0);

disp(datestr(now, 'dd/mm/yy-HH:MM'))

fprintf('E: %d (should be 260), Eafter: %d (should be 44)\n', E, Eafter);
fprintf('unique(labels) should be [0 4] and is: [');
fprintf('%d ', unique(labels));
fprintf(']\n')

disparityLabels = (reshape(labels,[param.W,param.H]));
disparityLabels = disparityLabels';
figure(100 + refImgIndex); subplot(1,2,2); imshow(rescale(disparityLabels)); title('Disparity Init')