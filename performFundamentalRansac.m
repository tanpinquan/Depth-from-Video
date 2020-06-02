function [outputF, inlierIdx]= performFundamentalRansac(xPoints,xPrimePoints,matches,nPoints,threshold,nRuns,bundleAdjust)
xMatchedPoints = xPoints(1:2,matches(1,:));
xPrimeMatchedPoints = xPrimePoints(1:2,matches(2,:));
xMatchedPoints(3,:) = 1;
xPrimeMatchedPoints(3,:) = 1;

points{1} = [107 306; 166 321; 444 329; 219 398; 571 373; 725 483; 782 388;151 514; 665 386];
points{2} = [216 296; 277 310; 582 326; 343 391; 675 372; 939 485; 917 391;315 508; 728 388];
points{1}(:,3) = 1;
points{2}(:,3) = 1;

maxInliers = 0;
outputF = eye(3);
inlierIdx = [];
for i =1:nRuns
    perm = randperm(size(matches,2)) ;
    selectedMatches = matches(:,perm(1:nPoints));
%     selectedMatches = matches(:,1:4);
    selectedXPoints = xPoints(1:2,selectedMatches(1,:))';
    selectedXPrimePoints = xPrimePoints(1:2,selectedMatches(2,:))';
    
    %% Plot
%     subplot(1,2,1);hold on;plot(selectedXPoints(:,1),selectedXPoints(:,2),'xb', 'markersize', 10)
%     subplot(1,2,2);hold on;plot(selectedXPrimePoints(:,1),selectedXPrimePoints(:,2),'xb', 'markersize', 10)

    F = computeFundamentalMat(selectedXPoints,selectedXPrimePoints);
    F2 = computeFundamentalMat(points{1},points{2});

    for j = 1:length(xMatchedPoints)
        errors(j) = abs(xPrimeMatchedPoints(:,j)'*F*xMatchedPoints(:,j));
%         disp(points{2}(i,:)*F{2}*points{1}(i,:)')
    end

    for j = 1:length(points{1})
        errors2(j) = abs(points{1}(j,:)*F*points{2}(j,:)');
    end
    
    numInliers = sum(errors<threshold);
%     inlierIdx = [];
    if(numInliers>maxInliers) 
%         i
        maxInliers = numInliers;
        display(['iteration ' num2str(i) ': ' num2str(maxInliers) ' inliers'])
        inlierIdx = find(errors<threshold);
        outputF = F;
        
    end
end

if (isempty(inlierIdx))
    outputF = [];

else
    xInliers = xMatchedPoints(:,inlierIdx)';
    xPrimeInliers = xPrimeMatchedPoints(:,inlierIdx)';

    
    outputF = computeFundamentalMat(xInliers,xPrimeInliers);
%     outputH = outputH./outputH(3,3);
    if(bundleAdjust)
        outputF = bundleAdjustmentFundamental(outputF,xPrimeInliers,xInliers);
    end
end




