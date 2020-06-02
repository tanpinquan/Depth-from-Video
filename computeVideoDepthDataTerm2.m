function [dataTerms] = computeVideoDepthDataTerm2(testLoc, mat1, mat2, refImgIndex, img, imgRGB,param)
% dataTerms = cell(1,nImages);
% diffTerms = cell(1,nImages);
% diffTerms = inf(nLabels,nImages);
startImage = max(1,refImgIndex - param.imageWindowSize);
endImage = min(param.nImages, refImgIndex + param.imageWindowSize);

nImages = endImage-startImage+1;
diffTerms = nan(param.nLabels,nImages);

x = cell(1,nImages);
i = refImgIndex;
for j = startImage:endImage
    if(i~=j)
        x{j} = mat1{i,j}*testLoc + param.d.*mat2{i,j};
        x{j} = round(x{j}./x{j}(3,:));
        validPoints = x{j}(1,:)>0 & x{j}(1,:)<=param.W & x{j}(2,:)>0 & x{j}(2,:)<=param.H;
        x{j} = x{j}(:,validPoints);

    %     testPixel2 = squeeze(img{i}(testLoc(2),testLoc(1),:));
        testPixel = permute(img{i}(testLoc(2),testLoc(1),:),[3 2 1]);

        idx = sub2ind([param.H param.W], x{j}(2,:), x{j}(1,:));
        pixelLine = [imgRGB.r{j}(idx); imgRGB.g{j}(idx); imgRGB.b{j}(idx)];


        %% Photometric
        colourDiff = (sum(abs(pixelLine - testPixel),1)/3)';

        if(length(colourDiff)>= param.nLabels)
            diffTerms(1:param.nLabels,j-startImage+1) = colourDiff(1:param.nLabels);
        else
            diffTerms(1:length(colourDiff),j-startImage+1) = colourDiff;
        end
        
        %% Debug
%         if(isequal(testLoc, [700; 260; 1]))
%             figure(700+j);
%             subplot(1,2,1); imshow(uint8(img{i})); title(['image ' num2str(i)]);
%             subplot(1,2,2); imshow(uint8(img{j})); title(['image ' num2str(j)]);
% 
%             subplot(1,2,1)
%             hold on;plot(testLoc(1),testLoc(2),'x')
% 
%             subplot(1,2,2)
%             hold on;plot(x{i,j}(1,:),x{i,j}(2,:),'x')          
%         end
        
    end


end

diffTerms = (param.sigma./(param.sigma + diffTerms));

% dataTerms = sum(diffTerms,2);

dataTerms = mean(diffTerms,2,'omitnan');
dataTerms(isnan(dataTerms))=0;

maxTerm = max(dataTerms);
if(maxTerm == 0)
    dataTerms(:) = 1;
else
    dataTerms = 1-dataTerms./max(dataTerms);
end

% if(isequal(testLoc, [700; 260; 1]))
%     disp('test'); 
%     [a idx] = max(diffTerms)
%     for j = 1:param.nImages
%         if(i~=j)
%             figure(700+j);
%             subplot(1,2,2)
%             hold on;plot(x{i,j}(1,idx(j)),x{i,j}(2,idx(j)),'rx')              
%         end
%     end
% end

    
