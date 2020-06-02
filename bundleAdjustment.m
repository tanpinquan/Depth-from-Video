function outputH = bundleAdjustment(h, initialPoints, targetPoints)

h=h./h(3,3);
initialPoints = initialPoints';
targetPoints = targetPoints';
h = h';


xData = initialPoints(1:2,:);
yData = targetPoints(1:2,:);

% test = applyHomography(h,xData);

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','off');
lb = [];
ub = [];
[outputH, resErr] = lsqcurvefit(@applyHomography,h,xData,yData,lb,ub,options);
outputH = outputH';

initErr = applyHomography(h,xData)-yData;
initErr = sum(sum(initErr.^2,2))
resErr
% test = applyHomography(hNew,xData);
% test2 = applyHomography(h,xData);
% errorsOld = sum((targetPoints(1:2,:)-test2).^2,2)
% errorsNew = sum((targetPoints(1:2,:)-test).^2,2)

% 
% X.x = initialPoints(1,:);
% X.y = initialPoints(2,:);
% 
% XPrime.x = targetPoints(1,:);
% XPrime.y = targetPoints(2,:);
% 
% mappedPoints = (h(1:2,:)*initialPoints)./(h(3,:)*initialPoints);
% mappedPoints(3,:) = 1;
% errors = sum((targetPoints-mappedPoints).^2,2)
% 
% XMap.x = mappedPoints(1,:);
% XMap.y = mappedPoints(2,:);
% 
% 
% % mappedPoints = h*initialPoints';
% % mappedPoints = mappedPoints./mappedPoints(3,:);
% % mappedPoints = mappedPoints';
% % errors = sum((targetPoints-mappedPoints).^2)
% 
% 
% 
% 
% numerator.one = h(1,:)*initialPoints;
% denominator.one = h(3,:)*initialPoints;
% 
% numerator.two = h(2,:)*initialPoints;
% denominator.two = h(3,:)*initialPoints;
% J = zeros(3,3);
% 
% 
% J(1,1) = sum((-2*X.x).*(XPrime.x - XMap.x).*(1./denominator.one));
% J(1,2) = sum((-2*X.y).*(XPrime.x - XMap.x).*(1./denominator.one));
% J(1,3) = sum((-2).*(XPrime.x - XMap.x).*(1./denominator.one));
% 
% J(2,1) = sum((-2*X.x).*(XPrime.y - XMap.y).*(1./denominator.two));
% J(2,2) = sum((-2*X.y).*(XPrime.y - XMap.y).*(1./denominator.two));
% J(2,3) = sum((-2).*(XPrime.y - XMap.y).*(1./denominator.two));
% 
% J(3,1) = sum((2*X.x).*(XPrime.x - XMap.x).*XMap.x./denominator.one ...
%        + (2*X.x).*(XPrime.y - XMap.y).*XMap.y./denominator.two);
% 
% J(3,2) = sum((2*X.y).*(XPrime.x - XMap.x).*XMap.x./denominator.one ...
%        + (2*X.y).*(XPrime.y - XMap.y).*XMap.y./denominator.two);
%    
% J(3,3) = sum((2).*(XPrime.x - XMap.x).*XMap.x./denominator.one ...
%        + (2).*(XPrime.y - XMap.y).*XMap.y./denominator.two); 
% 
% % J(1,1) = sum((-2*X.x).*(XPrime.x - (numerator.one./denominator.one)).*(1./denominator.one));
% % J(1,2) = sum((-2*X.y).*(XPrime.x - (numerator.one./denominator.one)).*(1./denominator.one));
% % J(1,3) = sum((-2).*(XPrime.x - (numerator.one./denominator.one)).*(1./denominator.one));
% % 
% % J(2,1) = sum((-2*X.x).*(XPrime.y - (numerator.two./denominator.two)).*(1./denominator.two));
% % J(2,2) = sum((-2*X.y).*(XPrime.y - (numerator.two./denominator.two)).*(1./denominator.two));
% % J(2,3) = sum((-2).*(XPrime.y - (numerator.two./denominator.two)).*(1./denominator.two));
% 
% % J(3,2) = sum((2*X.x).*(XPrime.x - (numerator.one./denominator.one)).*(numerator.one./(denominator.one.^2)) ...
% %        + (2*X.x).*(XPrime.y - (numerator.two./denominator.two)).*(numerator.two./(denominator.two.^2))); 
% % % 
% % J(3,2) = sum((2*X.y).*(XPrime.x - (numerator.one./denominator.one)).*(numerator.one./(denominator.one.^2)) ...
% %        + (2*X.y).*(XPrime.y - (numerator.two./denominator.two)).*(numerator.two./(denominator.two.^2))); 
% %    
% % J(3,3) = sum((2).*(XPrime.x - (numerator.one./denominator.one)).*(numerator.one./(denominator.one.^2)) ...
% %        + (2).*(XPrime.y - (numerator.two./denominator.two)).*(numerator.two./(denominator.two.^2)));  
%    
%    
% hNew = h -0.00001*J
% 
% mappedPointsNew = (hNew(1:2,:)*initialPoints)./(hNew(3,:)*initialPoints);
% mappedPointsNew(3,:) = 1;
% errorsNew = sum((targetPoints-mappedPointsNew).^2,2)
% 
% % mappedPoints = h'*initialPoints';
% % mappedPoints = mappedPoints./mappedPoints(3,:);
% % mappedPoints = mappedPoints';
% % errors = sum((targetPoints-mappedPoints).^2)
% % 
% % J(1,1) = -2*sum((targetPoints(:,1)-mappedPoints(:,1)).*initialPoints(:,1));
% % J(1,2) = -2*sum((targetPoints(:,1)-mappedPoints(:,1)).*initialPoints(:,2));
% % J(1,3) = -2*sum((targetPoints(:,1)-mappedPoints(:,1)));
% % 
% % J(2,1) = -2*sum((targetPoints(:,2)-mappedPoints(:,2)).*initialPoints(:,1))
% % J(2,2) = -2*sum((targetPoints(:,2)-mappedPoints(:,2)).*initialPoints(:,2))
% % J(2,3) = -2*sum((targetPoints(:,2)-mappedPoints(:,2)))
% % 
% % J(3,1) = -2*sum((targetPoints(:,1)-mappedPoints(:,1)).*initialPoints(:,1))
% % J(3,2) = -2*sum((targetPoints(:,1)-mappedPoints(:,1)).*initialPoints(:,2))
% % J(3,3) = -2*sum((targetPoints(:,1)-mappedPoints(:,1)))



function F = applyHomography(x,xdata)
    xdata(3,:) = 1;
    F = x*xdata;
    F = F./F(3,:);
    F = F(1:2,:);
% end

