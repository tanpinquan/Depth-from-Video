function projectedFeatures = applyHomographyToFeatures(H,features)

projectedFeatures = features(1:2,:);
projectedFeatures(3,:) = 1;

projectedFeatures = H'*projectedFeatures;
projectedFeatures = projectedFeatures./projectedFeatures(3,:);

if(size(features,1)>2)
    projectedFeatures(3:4,:) = features(3:4,:);
end
