function cleanH5file(row,col,field,SourceF)
first_tp   = 1;
H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
cellpath_name = ['/field' num2str(field) '/cellpath'];
sisterList_name = ['/field' num2str(field) '/sisterList'];
bg_name = ['/field' num2str(field) '/bg'];

fileattrib(fullfile(SourceF,H5filename),'+w');


% Load Original seed points

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
    H5F.close(fid);
    cellpathinfo = h5info(fullfile(SourceF,H5filename), cellpath_name);
    
    cellpath_mat = h5read(fullfile(SourceF,H5filename),cellpath_name,[1 1 1], [cellpathinfo.Dataspace.Size(1) cellpathinfo.Dataspace.Size(2) cellpathinfo.Dataspace.Size(3)]);
    
    for tp=first_tp
        cellpath{tp} = cellpath_mat(:,:,tp);
    end
end

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
    H5F.close(fid);
    sisterListinfo = h5info(fullfile(SourceF,H5filename), sisterList_name);
    sisterList_mat = h5read(fullfile(SourceF,H5filename),sisterList_name,[1 1 1], [sisterListinfo.Dataspace.Size(1) sisterListinfo.Dataspace.Size(2) sisterListinfo.Dataspace.Size(3)]);
    
    for tp=first_tp
        sisterList{tp} = sisterList_mat(:,:,tp);
    end
end

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,bg_name,'H5P_DEFAULT')
    H5F.close(fid);
    bginfo = h5info(fullfile(SourceF,H5filename), bg_name);
    bg_mat = h5read(fullfile(SourceF,H5filename),bg_name,[1 1 1], [bginfo.Dataspace.Size(1) bginfo.Dataspace.Size(2) bginfo.Dataspace.Size(3)]);
    
    for tp=first_tp
        bg{tp} = bg_mat(:,:,tp);
    end
end

% Reassigning data

c_tp=1;
cellpath_mat = -1*(ones(size(cellpath{c_tp},1),2,length(cellpath)));
sisterList_mat = -1*(ones(size(sisterList{c_tp},1),size(sisterList{c_tp},2),length(sisterList)));
bg_mat = -1*(ones(size(bg{c_tp},1),2,length(bg)));

for tp=first_tp
        cellpath_mat(:,:,tp) = cellpath{tp};
        sisterList_mat(:,:,tp) = sisterList{tp};
        bg_mat(:,:,tp) = bg{tp};
end

if exist(fullfile(SourceF,H5filename),'file')
    fileattrib(fullfile(SourceF,H5filename),'+w');
    
    fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
    if ~H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
        H5F.close(fid);
        display(['Initializing ' H5filename ':' cellpath_name]);
    else
        H5L.delete(fid,cellpath_name,'H5P_DEFAULT');
        display(['Overwriting ' H5filename ':' cellpath_name]);
        H5F.close(fid);
    end
end

h5create(fullfile(SourceF,H5filename), cellpath_name, [size(cellpath_mat,1), size(cellpath_mat,2), size(cellpath_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(cellpath_mat,2), size(cellpath_mat,3)], 'Deflate', 9);
h5write(fullfile(SourceF,H5filename), cellpath_name, cellpath_mat, [1 1 1], [size(cellpath_mat,1) size(cellpath_mat,2) size(cellpath_mat,3)]);

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' sisterList_name]);
else
    H5L.delete(fid,sisterList_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' sisterList_name]);
    H5F.close(fid);
end

h5create(fullfile(SourceF,H5filename), sisterList_name, [size(sisterList_mat,1), size(sisterList_mat,2), size(sisterList_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(sisterList_mat,2), size(sisterList_mat,3)], 'Deflate', 9);
h5write(fullfile(SourceF,H5filename), sisterList_name, sisterList_mat, [1 1 1], [size(sisterList_mat,1) size(sisterList_mat,2) size(sisterList_mat,3)]);

fid = H5F.open(fullfile(SourceF,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
if ~H5L.exists(fid,bg_name,'H5P_DEFAULT')
    H5F.close(fid);
    display(['Initializing ' H5filename ':' bg_name]);
else
    H5L.delete(fid,bg_name,'H5P_DEFAULT');
    display(['Overwriting ' H5filename ':' bg_name]);
    H5F.close(fid);
end

h5create(fullfile(SourceF,H5filename), bg_name, [size(bg_mat,1), size(bg_mat,2), size(bg_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(bg_mat,2), size(bg_mat,3)], 'Deflate', 9);
h5write(fullfile(SourceF,H5filename), bg_name, bg_mat, [1 1 1], [size(bg_mat,1) size(bg_mat,2) size(bg_mat,3)]);
disp(['Your data is saved to ' H5filename]);



