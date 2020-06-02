function outputH = bundleAdjustmentFundamental(F, x, xPrime)



xData = [x xPrime];
yData = zeros(length(x),1);

initErr = applyFundamental(F,xData)-yData;

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','off');
lb = [];
ub = [];
[outputH, resErr] = lsqcurvefit(@applyFundamental,F,xData,yData,lb,ub,options);
outputH = outputH';

initErr = applyFundamental(F,xData)-yData;
initErr = sum(sum(initErr.^2,2))
resErr


function F = applyFundamental(x,xdata)
    for i = 1:length(xdata)       
        F(i,1) = xdata(i,1:3)*x*xdata(i,4:6)';
    end
% end

