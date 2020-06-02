clc
clear
close all


nImages = 2;
param.imageWindowSize = 1;
imageStep = 20;
startImage = 0;
cameraTxt = readmatrix('images/cameras.txt');



K = [1221.2270770	0.0000000	479.5000000	
0.0000000	1221.2270770	269.5000000	
0.0000000	0.0000000	1.0000000];
%% Load Images
for i=1:nImages
   disp(['loading image ' num2str(startImage+(i-1)*imageStep)]);
%    disp((startImage+(i-1)*imageStep)*7+1);

   img{i} = double(imread(['images/test' num2str(startImage+(i-1)*imageStep, "%.4d") '.jpg'])); 
   imgRGB.r{i} = img{i}(:,:,1);
   imgRGB.g{i} = img{i}(:,:,2);
   imgRGB.b{i} = img{i}(:,:,3);

   camera{i} = cameraTxt((startImage+(i-1)*imageStep)*7+1:(startImage+(i-1)*imageStep)*7+7,1:3);
   K_ref{i} = camera{i}(1:3,:);
   R_ref{i} = camera{i}(4:6,:);
   T_ref{i} = camera{i}(7,:);   
     
    figure(1);
    subplot(1,nImages,i); imshow(uint8(img{i})) ;    

end

points{1} = [107 306; 166 321; 444 329; 219 398; 571 373; 725 483; 782 388;151 514; 665 386];
points{2} = [216 296; 277 310; 582 326; 343 391; 675 372; 939 485; 917 391;315 508; 728 388];


points{1} = [166 306; 445 328; 569 373; 724 482; 782 389; 152 518; 122 250; 57 460]
points{2} = [276 296; 582 326; 674 370; 939 484; 919 392; 308 509; 234 240; 151 449]

%% Fundamental matrix
[F{1}] = estimateFundamentalMatrix(points{1},points{2},'Method','Norm8Point');

[F{2}] = estimateFundamentalMatrix(points{2},points{1},'Method','Norm8Point');
% [F{2}] = computeFundamentalMat(points{1},points{2});

%% Epipolar line 1
epiLines = epipolarLine(F{1},points{1});
test1 = [163;318;1]

l = F{1}*test1;
for i=1:nImages
    figure(2);
    subplot(1,nImages,i)
    imshow(uint8(img{i}))
end

xStart = 1;
xStep = 1;

xTest = xStart:xStep:1000*xStep;

yTest = (-l(1).*xTest - l(3))/l(2);

subplot(1,2,1);hold on;plot(test1(1),test1(2),'yx', 'markersize',10','linewidth',2)
subplot(1,2,2);hold on;plot(xTest,yTest,'y.', 'markersize',10','linewidth',2)

%% Epipolar line 2
epiLines = epipolarLine(F{2},points{1});
test1 = [673;372;1]

l = F{2}*test1;


xStart = 1;
xStep = 1;

xTest = xStart:xStep:1000*xStep;

yTest = (-l(1).*xTest - l(3))/l(2);

subplot(1,2,2);hold on;plot(test1(1),test1(2),'rx', 'markersize',10','linewidth',2)
subplot(1,2,1);hold on;plot(xTest,yTest,'r.', 'markersize',10','linewidth',2)


%% Essesntial Matrix
% [E_est{1},R_est1{1},R_est2{1}, t_est1{1}, t_est2{1}] = computeEssentialMat(F{1},K,K);
[E{2},R2_est{1},R2_est{2}, T2_est{1}, T2_est{2}] = computeEssentialMat(F{2},K,K);

cameraParams = cameraParameters('IntrinsicMatrix',K); 

[relativeOrientation2,relativeLocation2] = relativeCameraPose(E{2},cameraParams,points{2},points{1})
[rotationMatrix2,translationVector2] = cameraPoseToExtrinsics(relativeOrientation2,relativeLocation2)

k = 1;
for i = 1:2
    for j= 1:2
        R2{k} = R2_est{i}
        T2{k} = T2_est{j}
        C2{k} = -R2{k}'*T2{k};
        k = k+1;
    end
end
% R2 = R2_est{2};
% T2 = t_est2{2};
% C2 = -R2'*T2;


% C2 = T_ref{2}';


R1 = [1 0 0; 0 1 0; 0 0 1];
T1 = [0 0 0]';
C1 = -R1'*T1;

%% Get epipolar line

mat1_ref = K*(R_ref{2}')*(R_ref{1})*inv(K);
mat2_ref = K*(R_ref{2}')*(T_ref{1}'-T_ref{2}');
mat1_ref*test1

%% Get 4 epipolar lines

for i=1:nImages
    figure(100);
    subplot(1,nImages,i)
    imshow(uint8(img{i}))
end

% test1 = [154;520;1];
% test1 = [166;306;1];
test1 = [725;485;1];
% test1 = [140;241;1];

plotMarker = ['r' 'b' 'y' 'k']
for i = 1:4
    mat1 = K*R2{i}*(R1')*inv(K);
    mat2 = K*R2{i}*(C1-C2{i});
    mat1*test1

    nLabels = 1000;
    d = 0:0.0002:nLabels*0.2;
    d = d(1:nLabels);




    x = mat1*test1 + d.*mat2;
    x = x./x(3,:);
    subplot(1,2,1);hold on;plot(test1(1),test1(2),'rx', 'markersize',10','linewidth',2)

    subplot(1,nImages,2)
    hold on;plot(x(1,:),x(2,:),'.', 'color', plotMarker(i), 'markersize',5','linewidth',3)
end

%% Ref
% x_ref = mat1_ref*test1 + d.*mat2_ref;
% x_ref = x_ref./x_ref(3,:);
% subplot(1,2,1);hold on;plot(test1(1),test1(2),'rx')
% 
% subplot(1,nImages,2)
% hold on;plot(x_ref(1,:),x_ref(2,:),'bx')
