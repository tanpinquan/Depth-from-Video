function [dataTerms] = computeUnrectifiedDepthDataTerm(x, d, mat1, mat2, img1, img2, W, H, nLabels, method, sigma)


xPrimes = mat1*x + d.*mat2;
xPrimes = xPrimes./xPrimes(3,:);

xPrimes = round(xPrimes);

validPoints = xPrimes(1,:)>0 & xPrimes(1,:)<=W & xPrimes(2,:)>0 & xPrimes(2,:)<=H;
xPrimes = xPrimes(:,validPoints);

testPixel = squeeze(img1(x(2),x(1),:));



idx = sub2ind(size(img2.r), xPrimes(2,:), xPrimes(1,:));

pixelLine = [img2.r(idx); img2.g(idx); img2.b(idx)];

colourDiff = sum(abs(pixelLine - testPixel),1)/3;

if(strcmp(method,'norm'))
    colourDiff = sigma./(sigma + colourDiff(:,2:end));
    colourDiff = 1-colourDiff./max(colourDiff);
    dataTerms = ones(nLabels,1);

else
    dataTerms = 999*ones(nLabels,1);
end

if(length(colourDiff)>= nLabels)
    dataTerms(1:nLabels) = colourDiff(1:nLabels);
else
    dataTerms(1:length(colourDiff)) = colourDiff;
end

% dataTerms = dataTerms';



