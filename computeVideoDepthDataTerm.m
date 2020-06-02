function [dataTerms] = computeVideoDepthDataTerm(testLoc, d, mat1, mat2, refImgIndex, img, imgRGB, W, H, nImages, nLabels, sigma)
% dataTerms = cell(1,nImages);
% diffTerms = cell(1,nImages);
% diffTerms = inf(nLabels,nImages);
diffTerms = nan(nLabels,nImages);

x = cell(1,nImages);
i = refImgIndex;
for j = 1:nImages
    if(i~=j)
        x{i,j} = mat1{i,j}*testLoc + d.*mat2{i,j};
        x{i,j} = round(x{i,j}./x{i,j}(3,:));
        validPoints = x{i,j}(1,:)>0 & x{i,j}(1,:)<=W & x{i,j}(2,:)>0 & x{i,j}(2,:)<=H;
        x{i,j} = x{i,j}(:,validPoints);

    %     testPixel2 = squeeze(img{i}(testLoc(2),testLoc(1),:));
        testPixel = permute(img{i}(testLoc(2),testLoc(1),:),[3 2 1]);

        idx = sub2ind([H W], x{i,j}(2,:), x{i,j}(1,:));
        pixelLine = [imgRGB.r{j}(idx); imgRGB.g{j}(idx); imgRGB.b{j}(idx)];


        %% Photometric
        colourDiff = (sum(abs(pixelLine - testPixel),1)/3)';

        if(length(colourDiff)>= nLabels)
            diffTerms(1:nLabels,j) = colourDiff(1:nLabels);
        else
            diffTerms(1:length(colourDiff),j) = colourDiff;
        end
        
    end


end
diffTerms = (sigma./(sigma + diffTerms));

% dataTerms = sum(diffTerms,2);

dataTerms = mean(diffTerms,2,'omitnan');
dataTerms(isnan(dataTerms))=0;

maxTerm = max(dataTerms);
if(maxTerm == 0)
    dataTerms(:) = 1;
else
    dataTerms = 1-dataTerms./max(dataTerms);
end
    
    
