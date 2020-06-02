function P = computeProjectiveMat(F)

[U,D,V] = svd(F);
Z = [0 1 0; -1 0 0; 0 0 0];

tSkew = U*Z*U';
M = U*Z'*D*V';

t = [tSkew(3,2); tSkew(1,3); tSkew(2,1)]

P = [t M];

