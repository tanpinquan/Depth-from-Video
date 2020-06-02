function [E,R1,R2,t1,t2] = computeEssentialMat(F,K,KPrime)

E = KPrime'*F*K;

[U,D,V] = svd(E);
Z = [0 1 0; -1 0 0; 0 0 0];
W = [0 -1 0; 1 0 0; 0 0 1];
Winv = W';
% 
tSkew = U*Z*U';
tSkew2 = U*W*D*U';

R1 = U*Winv*V';
R2 = U*W*V';
R1 = det(R1)*R1;
R2 = det(R2)*R2;


t1 = [tSkew(3,2); tSkew(1,3); tSkew(2,1)];
% t1 = [D(1,1) D(2,2) 1]'.*t1;
t2 = -t1;

t1 = [tSkew2(3,2); tSkew2(1,3); tSkew2(2,1)];
t2 = -t1;
