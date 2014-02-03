function [bwNuc, bwOther] = detectSeededNuclei(im, nucCoord, avgNucDiameter)

minNucArea = round(0.1*avgNucDiameter^2);
maxMAL = 5*avgNucDiameter;
minOtherNucArea = round(0.25*avgNucDiameter^2);
padWidth = 10;

%% pad image to facilitate detection of objects on the image boundary
[ySize, xSize] = size(im);
im=im2uint16(im);
[iMax, iWidth] = contrastValues(im);
tmp = iMax + (rand(ySize+2*padWidth, xSize+2*padWidth)-0.5) * 2*iWidth;
tmp(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth) = im;
im = tmp/max(tmp(:));

%% convert nuclear coordinates to a binary mask
nucCent = false(size(im));
for n = 1: size(nucCoord, 1)
    if nucCoord(n,1) > 0 && nucCoord(n,2) > 0
        nucCent(nucCoord(n,2)+padWidth, nucCoord(n,1)+padWidth) = true;
    end
end

%% edge detection
% convert image to double
im = im2double(im);
% smooth image
imS = filterGauss2D(im,1);
% call steerable detector to detect nuclear boundary
[~, ~, nms, ~] = steerableDetector(im, 3, 2);
% normalize the non-maximum suppressed results by the original image
nmsNorm = nms ./ im;

%% reconstruct nuclei from edges
bwNuc = false(size(im));
bwOther = false(size(im));
% loop over different level of threshold of the NMS
quantileCutoff = 0.99;
while quantileCutoff > 0.96
    % threshold normalize NMS
    bwEdgeRaw = nmsNorm > quantile(nmsNorm(:), quantileCutoff);
    % loop over different dilation radius to reconstruct nuclei
    r = 0;
    while r < 4
        if r > 0
            bwEdge = imdilate(bwEdgeRaw, strel('disk', r));
            bwNucTmp = imfill(bwEdge, 'holes');
            bwNucTmp = imerode(bwNucTmp, strel('disk', r));
        else
            bwEdge = bwEdgeRaw;
            bwNucTmp = imfill(bwEdge, 'holes');
        end
        % remove random bits on the peripheral of nuclei
        bwNucTmp = imopen(bwNucTmp, strel('disk', 2));
        % ensure new objects not touching the existing nuclear objects
        tmp=bwdist(bwNuc);
        tmp=tmp<2;
        bwNucTmp = ~tmp & bwNucTmp;
        
        % divide up objects that contain more than one nuclear center
% %         distMat = -bwdist(~bwNucTmp);
%         distMat = -imS;
%         distMat = imimposemin(distMat, nucCent);
%         lNucTmp = watershed(distMat);
%         bwNucTmp = bwNucTmp & (lNucTmp > 0);        
        [lNucTmp, nNucTmp] = bwlabel(bwNucTmp);
        bwNucTmp = false(size(bwNuc));
        for n = 1:nNucTmp
            nCent = sum(sum(nucCent(lNucTmp == n)));
            if nCent <= 1
                bwNucTmp = bwNucTmp | lNucTmp == n;
            elseif nCent > 1
                bw = lNucTmp == n;
                distMat = bwdist(nucCent&bw);
                distMat = imimposemin(distMat, nucCent);
                l = watershed(distMat);
                bw = bw & (l > 0);
                bwNucTmp = bwNucTmp | bw;
            end
        end
                
        
        % filter out non-nuclear objects
        if max(bwNucTmp(:)) > 0
            ffactorCutoff = 0.9 - 0.1*r;
            [lNucTmp, nNucTmp] = bwlabel(bwNucTmp);
            S = regionprops(lNucTmp,'Perimeter','Area','MajorAxisLength');
            for n = 1: nNucTmp
                if S(n).MajorAxisLength > maxMAL || S(n).Area < minNucArea || 4*pi*S(n).Area/S(n).Perimeter^2 < ffactorCutoff
                    bwNucTmp(lNucTmp == n) = false;
                end
            end
        end
        
        % keep objects with corresponding nuclear centers
        lNucTmp = bwlabel(bwNucTmp);
        idx = unique(lNucTmp(nucCent));
        idx(idx == 0) = [];
        bwNucTmp = ismember(lNucTmp, idx);
        bwOtherTmp = lNucTmp > 0 & ~bwNucTmp;
        
        % update detected nuclear mask
        bwNuc = bwNuc | bwNucTmp;
        if r <= 1
            bwOther = bwOther | bwOtherTmp;
        end
        % update nuclear centers
        nucCent = nucCent & (~bwNucTmp);
        
        if max(nucCent(:)) > 0
            r = r + 1;
        else
            r = 4;
        end
    end
    if max(nucCent(:)) > 0
        quantileCutoff = quantileCutoff - 0.01;
    else
        quantileCutoff = 0.96;
    end
end

% recover seeded nuclei that was not detected by edge
if max(nucCent(:)) > 0
    bw = blobSegmentThresholdGeneral(im, 'rosin', 1, 1, minNucArea);
    tmp = bwdist(bwNuc);
    tmp = tmp < 2;
    bw = ~tmp & bw;
    l = bwlabel(bw);
    idx = unique(l(nucCent));
    idx(idx == 0) = [];
    bwNucTmp = ismember(l, idx);
    
    if max(bwNucTmp(:)) > 0
        [lNucTmp, nNucTmp] = bwlabel(bwNucTmp);
        bwNucTmp = false(size(bwNuc));
        for n = 1:nNucTmp
            nCent = sum(sum(nucCent(lNucTmp == n)));
            if nCent <= 1
                bwNucTmp = bwNucTmp | lNucTmp == n;
            elseif nCent > 1
                bw = lNucTmp == n;
                distMat = bwdist(nucCent&bw);
                distMat = imimposemin(distMat, nucCent);
                l = watershed(distMat);
                bw = bw & (l > 0);
                bwNucTmp = bwNucTmp | bw;
            end
        end
        bwNuc = bwNuc | bwNucTmp;
        nucCent = nucCent & (~bwNuc);
    end
end

% % place a sphere in the spot of missing nuclei
% if max(nucCent(:)) > 0
%     bwNucTmp = imdilate(nucCent, strel('disk',round(avgNucDiameter/2)));
%     tmp=bwdist(bwNuc);
%     tmp=tmp<2;
%     bwNucTmp = ~tmp & bwNucTmp;
% end

% recover non-seeded nuclei
% bwOtherTmp = blobSegmentThresholdGeneral(im, 'rosin', 1, 1, minNucArea);
% lOtherTmp = bwlabel(bwOtherTmp);
% idx = unique(lOtherTmp(bwNuc | bwOther));
% idx(idx == 0) = [];
% tmp = ismember(lOtherTmp, idx);
% bwOtherTmp = bwOtherTmp & ~tmp;
% bwOther = bwOther | bwOtherTmp;
% bwOther = bwareaopen(bwOther, minOtherNucArea);

% recover non-seeded nuclei
lOther=bwlabel(bwOther);
idx=unique(lOther(bwNuc));
idx(idx==0)=[];
bwOther = bwOther & ~ismember(lOther,idx);
bwOther = bwareaopen(bwOther,minOtherNucArea);

% remove padded edge region
bwNuc = bwNuc(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth);
bwOther = bwOther(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth);

end