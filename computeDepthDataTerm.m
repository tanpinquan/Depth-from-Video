function [dataTerms] = computeDepthDataTerm(data1, pixelRow, index, maxDisplacement, method, sigma)

    if(index==101)
        index;
    end
    pixelRow = squeeze(pixelRow);
    data1 = squeeze(data1);
    rDiff = pixelRow-data1(1);
    colourDiff = pixelRow - data1';
    
    
    
    colourDiff = sum(abs(pixelRow - data1'),2)/3;
    
%     dataTerms(W-index:W-index+W-1) = colourDiff;
%     dataTerms = dataTerms';
    
    startIndex = max(1, index-maxDisplacement+1);
    endIndex = min(length(pixelRow), index+maxDisplacement+1);
    
    colourDiff = sum(abs(pixelRow(startIndex:endIndex,:) - data1'),2)/3;
    
%     dataTerms = 999*ones(1,2*maxDisplacement+1);
    if(strcmp(method,'norm'))
        colourDiff = sigma./(sigma + colourDiff);
        colourDiff = 1-colourDiff./max(colourDiff);
        dataTerms = ones(1,2*maxDisplacement+1,1);

    else
        dataTerms = 999*ones(1,2*maxDisplacement+1,1);
    end
    
    startDataTermIndex = maxDisplacement + startIndex-(index);
    endDataTermIndex = maxDisplacement + endIndex-(index);
    dataTerms(startDataTermIndex:endDataTermIndex) = colourDiff;
    dataTerms = dataTerms';
    
end
