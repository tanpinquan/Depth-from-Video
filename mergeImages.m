function [outputImages, canvasRef] = mergeImages(image1,ref1,image2,ref2,mode)


ref1.XWorldLimits = floor(ref1.XWorldLimits);
ref1.YWorldLimits = floor(ref1.YWorldLimits);
ref2.XWorldLimits = round(ref2.XWorldLimits);
ref2.YWorldLimits = round(ref2.YWorldLimits);

canvasRef = imref2d;
canvasRef.XWorldLimits(1) = min(ref1.XWorldLimits(1),ref2.XWorldLimits(1));
canvasRef.XWorldLimits(2) = max(ref1.XWorldLimits(2),ref2.XWorldLimits(2));
canvasRef.YWorldLimits(1) = min(ref1.YWorldLimits(1),ref2.YWorldLimits(1));
canvasRef.YWorldLimits(2) = max(ref1.YWorldLimits(2),ref2.YWorldLimits(2));
canvasOffset.x = -(min(canvasRef.XWorldLimits));
canvasOffset.y = -(min(canvasRef.YWorldLimits));
canvas1 = double(zeros(canvasRef.ImageExtentInWorldY,canvasRef.ImageExtentInWorldX,3));
canvas2 = double(zeros(canvasRef.ImageExtentInWorldY,canvasRef.ImageExtentInWorldX,3));

canvas1Limits.x = [1+ref1.XWorldLimits(1) + canvasOffset.x, ref1.XWorldLimits(2) + canvasOffset.x];
canvas1Limits.y = [1+ref1.YWorldLimits(1) + canvasOffset.y, ref1.YWorldLimits(2) + canvasOffset.y];

canvas2Limits.x = [1+ref2.XWorldLimits(1) + canvasOffset.x, ref2.XWorldLimits(2) + canvasOffset.x];
canvas2Limits.y = [1+ref2.YWorldLimits(1) + canvasOffset.y, ref2.YWorldLimits(2) + canvasOffset.y];

canvas1(canvas1Limits.y(1):canvas1Limits.y(2),canvas1Limits.x(1):canvas1Limits.x(2),:) = image1;
canvas2(canvas2Limits.y(1):canvas2Limits.y(2),canvas2Limits.x(1):canvas2Limits.x(2),:) = image2;

if(strcmp(mode,'average'))
    mask1 = canvas1(:,:,1)==0 & canvas1(:,:,2) ==0 & canvas1(:,:,3) ==0;
    mask1 = repmat(mask1,1,1,3);
    canvas1(mask1) = nan;

    mask2 = canvas2(:,:,1)==0 & canvas2(:,:,2) ==0 & canvas2(:,:,3) ==0;
    mask2 = repmat(mask2,1,1,3);
    canvas2(mask2) = nan;

    r = mean(cat(3,canvas1(:,:,1),canvas2(:,:,1)),3,'omitnan');
    g = mean(cat(3,canvas1(:,:,2),canvas2(:,:,2)),3,'omitnan');
    b = mean(cat(3,canvas1(:,:,3),canvas2(:,:,3)),3,'omitnan');
    outputImages(:,:,1) = r;
    outputImages(:,:,2) = g;
    outputImages(:,:,3) = b;
    outputImages = uint8(outputImages);
elseif(strcmp(mode,'image1'))
    mask = canvas1(:,:,1)~=0 & canvas1(:,:,2) ~=0 & canvas1(:,:,3) ~=0;
    mask = repmat(mask,1,1,3);
    canvas2(mask) = canvas1(mask);
    outputImages = uint8(canvas2);
elseif(strcmp(mode,'image2'))
    mask = canvas2(:,:,1)~=0 & canvas2(:,:,2) ~=0 & canvas2(:,:,3) ~=0;
    mask = repmat(mask,1,1,3);
    canvas1(mask) = canvas2(mask);
    outputImages = uint8(canvas1); 
end

