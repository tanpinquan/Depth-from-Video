function [dataTerms] = computeBundleAdjustDataTerm2(testLoc, mat1, mat2, disparityInit, refImgIndex, img, imgRGB, param)


startImage = max(1,refImgIndex - param.imageWindowSize);
endImage = min(param.nImages, refImgIndex + param.imageWindowSize);

nImages = 2*param.imageWindowSize+1;
diffTerms = nan(param.nLabels,nImages);
distTerms = nan(param.nLabels,nImages);

x = cell(1,nImages);
xBack = cell(1,nImages);

i = refImgIndex;
for j = startImage:endImage
    if(i~=j)
        
        x{j} = mat1{i,j}*testLoc + param.d.*mat2{i,j};
        x{j} = round(x{j}./x{j}(3,:));
        validPoints = x{j}(1,:)>0 & x{j}(1,:)<=param.W & x{j}(2,:)>0 & x{j}(2,:)<=param.H;
        x{j} = x{j}(:,validPoints);

        %% Photometric
        testPixel = permute(img{i}(testLoc(2),testLoc(1),:),[3 2 1]);

        idx = sub2ind([param.H param.W], x{j}(2,:), x{j}(1,:));
        pixelLine = [imgRGB.r{j}(idx); imgRGB.g{j}(idx); imgRGB.b{j}(idx)];

        colourDiff = (sum(abs(pixelLine - testPixel),1)/3)';

        
        %% Geometric
        disparityInd = disparityInit{j}(idx)+1;
        
        xBack{j} = mat1{j,i}*x{j} + param.d(disparityInd).*mat2{j,i};
        xBack{j} = xBack{j}./xBack{j}(3,:);

        distDiff = sum((xBack{j} - testLoc).^2);
        
        %% Store terms
        nValidTerms = length(colourDiff);
        if(nValidTerms>= param.nLabels)
            diffTerms(1:param.nLabels,j-startImage+1) = colourDiff(1:param.nLabels);
            distTerms(1:param.nLabels,j-startImage+1) = distDiff(1:param.nLabels);
        else
            diffTerms(1:length(colourDiff),j-startImage+1) = colourDiff;
            distTerms(1:length(colourDiff),j-startImage+1) = distDiff;
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
    
    
