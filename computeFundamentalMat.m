function F = computeFundamentalMat(x,xPrime)
% x = x';
% xPrime = xPrime';
nPoints = size(x,1);

% for i = 1:nPoints
%     A(i,:) = [xPrime(i,1)*x(i,1) xPrime(i,1)*x(i,2) xPrime(i,1) xPrime(i,2)*x(i,1) xPrime(i,2)*x(i,2) xPrime(i,2) x(i,1) x(i,2) 1];
%     
% %     A2((i-1)*2+1,:) = [initialPoints(i,1) initialPoints(i,2) 1 0 0 0 -mappedPoints(i,1)*initialPoints(i,1) -mappedPoints(i,1)*initialPoints(i,2)];
% %     A2((i-1)*2+2,:) = [0 0 0 initialPoints(i,1) initialPoints(i,2) 1 -mappedPoints(i,2)*initialPoints(i,1) -mappedPoints(i,2)*initialPoints(i,2)]; 
% end

A = [xPrime(:,1).*x(:,1) xPrime(:,1).*x(:,2) xPrime(:,1) xPrime(:,2).*x(:,1) xPrime(:,2).*x(:,2) xPrime(:,2) x(:,1) x(:,2)];
A(:,9) = 1


[U,D,V] = svd(A);
F = V(:,end);
F = reshape(F,[3 3]);



% b = reshape(mappedPoints',[],1);

% h2 = inv(A2'*A2)*A2'*b;
% h2(end+1) = 1;
% h2 = reshape(h2,[3 3]);



