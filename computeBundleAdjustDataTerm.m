function [dataTerms] = computeBundleAdjustDataTerm(testLoc, mat1, mat2, disparityInit, refImgIndex, img, imgRGB, param)
% dataTerms = cell(1,nImages);
% diffTerms = cell(1,nImages);
% diffTerms = inf(nLabels,nImages);
diffTerms = nan(param.nLabels,param.nImages);
distTerms = nan(param.nLabels,param.nImages);

x = cell(param.nImages,param.nImages);
xBack = cell(param.nImages,param.nImages);

i = refImgIndex;
for j = 1:param.nImages
    if(i~=j)
        
        x{i,j} = mat1{i,j}*testLoc + param.d.*mat2{i,j};
        x{i,j} = round(x{i,j}./x{i,j}(3,:));
        validPoints = x{i,j}(1,:)>0 & x{i,j}(1,:)<=param.W & x{i,j}(2,:)>0 & x{i,j}(2,:)<=param.H;
        x{i,j} = x{i,j}(:,validPoints);

        %% Photometric
        testPixel = permute(img{i}(testLoc(2),testLoc(1),:),[3 2 1]);

        idx = sub2ind([param.H param.W], x{i,j}(2,:), x{i,j}(1,:));
        pixelLine = [imgRGB.r{j}(idx); imgRGB.g{j}(idx); imgRGB.b{j}(idx)];

        colourDiff = (sum(abs(pixelLine - testPixel),1)/3)';

        
        %% Geometric
        disparityInd = disparityInit{j}(idx)+1;
        
        xBack{i,j} = mat1{j,i}*x{i,j} + param.d(disparityInd).*mat2{j,i};
        xBack{i,j} = xBack{i,j}./xBack{i,j}(3,:);

        distDiff = sum((xBack{i,j} - testLoc).^2);
        
        %% Store terms
        nValidTerms = length(colourDiff);
        if(nValidTerms>= param.nLabels)
            diffTerms(1:param.nLabels,j) = colourDiff(1:param.nLabels);
            distTerms(1:param.nLabels,j) = distDiff(1:param.nLabels);
        else
            diffTerms(1:length(colourDiff),j) = colourDiff;
            distTerms(1:length(colourDiff),j) = distDiff;
        end
        
%         figure();
%         subplot(1,2,1); imshow(uint8(img{i})); title(['image ' num2str(i)]);
%         subplot(1,2,2); imshow(uint8(img{j})); title(['image ' num2str(j)]);
% 
%         subplot(1,2,1)
%         hold on;plot(testLoc(1),testLoc(2),'x')
% 
%         subplot(1,2,2)
%         hold on;plot(x{i,j}(1,:),x{i,j}(2,:),'x')       
%         
%         subplot(1,2,1)
%         hold on;plot(xBack{i,j}(1,:),xBack{i,j}(2,:),'x')       
    end


end
% figure(210+refImgIndex+1);surf(flipud(disparityInit{i}),'EdgeColor','None');
% colormap('gray');
% view(2);

diffTerms = (param.sigma./(param.sigma + diffTerms));
distTerms = exp(-distTerms./(2*param.sigmaDistSq));
% dataTerms = sum(diffTerms,2);
combinedTerms = diffTerms.*distTerms;
dataTerms = mean(combinedTerms,2,'omitnan');
dataTerms(isnan(dataTerms))=0;

maxTerm = max(dataTerms);
if(maxTerm == 0)
    dataTerms(:) = 1;
else
    dataTerms = 1-dataTerms./max(dataTerms);
end
    
    
