function [matches,score] = matchDescriptors(d1,d2,thresh)
d1 = int16(d1);
d2 = int16(d2);
n1 = size(d1,2);
n2 = size(d2,2);
% matches = zeros(2,n1);
% matches(1,:) = 1:n1;
matches = zeros(2,0);
score = [];
for i=1:n1
    dist = d2-d1(:,i);
    dist = sum(dist.^2,1);
%     [score(i), matches(2,i)] = min(dist);    
    


    [minDist, minIdx] = min(dist);    
    if(sum(minDist*thresh>dist)<=1)
        matches(1,end+1) = i;
        matches(2,end) = minIdx;
        score(end+1) = minDist;
    end

end

[score, sortInd] = sort(score,'ascend');
matches = matches(:,sortInd);