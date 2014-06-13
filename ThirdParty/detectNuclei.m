function [combined_bw] = detectNuclei(im,minNucDiameter,maxNucDiameter)
padWidth =10;
[ySize, xSize] = size(im);
im=im2uint16(im);
[iMax, iWidth] = contrastValues(im);
tmp = iMax + (rand(ySize+2*padWidth, xSize+2*padWidth)-0.5) * 2*iWidth;
tmp(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth) = im;
im = tmp/max(tmp(:));

minFormfactor=0.7;

maxNucArea=round(pi*maxNucDiameter^2/4);

combined_bw = false(size(tmp));
BW_GoodCells = [];
testSet = linspace(1.1,4,5);

for s = 1:length(testSet)
    
    [im2,~]  = edge(tmp,'canny',0,testSet(s));
    im3 = ~im2;
    bw=im3-bwareaopen(im3,maxNucArea,4);
    bw=imfill(bw,'holes');
    CC = bwconncomp(bw, 8);
    L = labelmatrix(CC);
    S = regionprops(CC,'EquivDiameter','Area','Perimeter');
    nucArea=cat(1,S.Area);
    nucPerim=cat(1,S.Perimeter);
    nucEquiDiameter=cat(1,S.EquivDiameter);
    nucFormfactor=4*pi*nucArea./(nucPerim.^2);
    
    GoodCells = find(nucEquiDiameter >= minNucDiameter & ...
        nucEquiDiameter < maxNucDiameter & ...
        nucFormfactor>minFormfactor);
    
    BW_GoodCells{s} = ismember(L,GoodCells);
    combined_bw = combined_bw | ismember(L,GoodCells);
end


% [L, num] = bwlabel(combined_bw, 8);
% newcombined_bw = false(size(combined_bw));
% for i=1:num
%     ObjectNo = zeros(length(testSet),1);
%     
%     for s = 1:length(testSet)
% 
%         CellBW = BW_GoodCells{s}.*ismember(L,i);
%         [~, num] = bwlabel(CellBW, 8);
%         ObjectNo(s) = num;
%     end
%     combined_cell = false(size(combined_bw));
%     if max(ObjectNo) > 1
%         GoodPlanes = find(ObjectNo>1,1);
%         combined_cell = combined_cell | (BW_GoodCells{GoodPlanes}.*ismember(L,i)  ) ;
%     elseif max(ObjectNo) == 1
%         GoodPlanes = find(ObjectNo==1,1);
%         combined_cell = combined_cell | (BW_GoodCells{GoodPlanes}.*ismember(L,i)  ) ;
%     end
%     newcombined_bw = newcombined_bw | combined_cell;
% end

% remove padded edge region
combined_bw = combined_bw(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth);

%showresult
im = im(padWidth+1:ySize+padWidth, padWidth+1:xSize+padWidth);
p2 = bwperim(combined_bw);
p2 = bwmorph(p2,'dilate',1);
im=imadjust(im2double(im));
imshow(cat(3,im,max(im-p2,p2),im ));

