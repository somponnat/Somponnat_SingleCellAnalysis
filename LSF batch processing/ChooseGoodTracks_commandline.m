function ChooseGoodTracks_commandline(SourceF,row,col,field,channel,fileformat,channelnames)

currentPath = pwd;
eval('cd ..');
addpath(genpath([pwd filesep 'ThirdParty']),'-end');
cd(currentPath);

nucFolder = 'nuclearMask';
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
fileattrib(fullfile(SourceF,H5filename),'+w');
cellpath_name = ['/field' num2str(field) '/cellpath'];
cellpathinfo = h5info(fullfile(SourceF,H5filename), cellpath_name);
fileattrib(fullfile(SourceF,H5filename),'+w');
fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
    H5F.close(fid);
    cellpath_mat = h5read(fullfile(SourceF,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
else
    return
end
totalTp = size(cellpath_mat,3);
presenceMat = zeros(size(cellpath_mat,1),totalTp);
for tp = 1:totalTp
    nuc_im   = imread(fullfile(SourceF,nucFolder, sprintf(fileformat,channelnames{channel},tp)));
    P = impixel(nuc_im,cellpath_mat(:,1,tp),cellpath_mat(:,2,tp));
    Pt = P(:,1);
    Pt(isnan(Pt))=0;
    presenceMat(:,tp) = Pt;
end

goodCells = [];
for cell = 1:size(cellpath_mat,1)
    if length(find(presenceMat(cell,:))) > 0.8*size(cellpath_mat,3)
        goodCells = [goodCells,cell];
    end
end

newCellpath_mat = cellpath_mat(goodCells,:,:);

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
H5L.delete(fid,cellpath_name,'H5P_DEFAULT');
display(['Overwriting ' H5filename ':' cellpath_name]);
H5F.close(fid);

h5create(fullfile(SourceF,H5filename), cellpath_name, ...
    [size(newCellpath_mat,1), size(newCellpath_mat,2), size(newCellpath_mat,3)], ...
    'Datatype', 'double', 'ChunkSize', [1, size(newCellpath_mat,2), size(newCellpath_mat,3)], 'Deflate', 9);
h5write(fullfile(SourceF,H5filename), cellpath_name, newCellpath_mat, [1 1 1], ...
    [size(newCellpath_mat,1) size(newCellpath_mat,2) size(newCellpath_mat,3)]);
