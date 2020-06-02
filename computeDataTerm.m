function [dataTerms] = computeDataTerm(data1, sourceColour, sinkColour)
    data1 = squeeze(data1);
    distSource = sum(abs(data1-sourceColour))/3;
    distSink = sum(abs(data1-sinkColour))/3;
    
%     if(distSink<45)
%         distSink = 45;
%     end
%      if(distSource<45)
%         distSource = 45;
%     end   
    dataTerms = [distSource;distSink];
    
    
%     dataTerms = round(dataTerms);
end
