clear all; close all;

rows = 4;
cols = 4:10;
field = 1;
plane = 1;
CHs= [1 4];
first_tp   = 1;
last_tp   = 315;
cellsize = 40;

for row=rows
    for col=cols
        for templateCH = CHs
            H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
            cellpath_name = ['/field' num2str(field) '/cellpath'];
            sisterList_name = ['/field' num2str(field) '/sisterList'];
            bg_name = ['/field' num2str(field) '/bg'];
            celltrackOUTfilename = ['celltrackOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];
            
            
            datasetname = ['/field' num2str(field) '/segmentsCH' num2str(templateCH)];
            selectedcells_name = ['/field' num2str(field) '/selectedcells'];
            maskOUTfilename = ['maskOUT_r' num2str(row) '_c' num2str(col) '_f' num2str(field) '_p' num2str(plane) '_ch' num2str(templateCH)];
            
            if exist([celltrackOUTfilename '.mat'],'file')
                
                load(celltrackOUTfilename);
                
                cellpath_mat = -1*(ones(size(cellpath{last_tp},1),2,last_tp));
                sisterList_mat = -1*(ones(size(sisterList{last_tp},1),size(sisterList{last_tp},2),last_tp));
                bg_mat = -1*(ones(size(bg{last_tp},1),2,last_tp));
                
                for tp=first_tp:last_tp
                    cellpath_mat(:,:,tp) = cellpath{tp};
                    sisterList_mat(:,:,tp) = sisterList{tp};
                    bg_mat(:,:,tp) = bg{tp};
                end
                
                if exist(H5filename,'file')
                    fileattrib(H5filename,'+w');
                    
                    fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                    if ~H5L.exists(fid,cellpath_name,'H5P_DEFAULT')
                        H5F.close(fid);
                        display(['Initializing ' H5filename ':' cellpath_name]);
                    else
                        H5L.delete(fid,cellpath_name,'H5P_DEFAULT');
                        display(['Overwriting ' H5filename ':' cellpath_name]);
                        H5F.close(fid);
                    end
                else
                    display(['Initializing ' H5filename ':' cellpath_name]);
                end
                
                h5create(H5filename, cellpath_name, [size(cellpath_mat,1), size(cellpath_mat,2), size(cellpath_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(cellpath_mat,2), size(cellpath_mat,3)], 'Deflate', 9);
                h5write(H5filename, cellpath_name, cellpath_mat, [1 1 1], [size(cellpath_mat,1) size(cellpath_mat,2) size(cellpath_mat,3)]);
                
                fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                if ~H5L.exists(fid,sisterList_name,'H5P_DEFAULT')
                    H5F.close(fid);
                    display(['Initializing ' H5filename ':' sisterList_name]);
                else
                    H5L.delete(fid,sisterList_name,'H5P_DEFAULT');
                    display(['Overwriting ' H5filename ':' sisterList_name]);
                    H5F.close(fid);
                end
                
                h5create(H5filename, sisterList_name, [size(sisterList_mat,1), size(sisterList_mat,2), size(sisterList_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(sisterList_mat,2), size(sisterList_mat,3)], 'Deflate', 9);
                h5write(H5filename, sisterList_name, sisterList_mat, [1 1 1], [size(sisterList_mat,1) size(sisterList_mat,2) size(sisterList_mat,3)]);
                
                fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                if ~H5L.exists(fid,bg_name,'H5P_DEFAULT')
                    H5F.close(fid);
                    display(['Initializing ' H5filename ':' bg_name]);
                else
                    H5L.delete(fid,bg_name,'H5P_DEFAULT');
                    display(['Overwriting ' H5filename ':' bg_name]);
                    H5F.close(fid);
                end
                
                h5create(H5filename, bg_name, [size(bg_mat,1), size(bg_mat,2), size(bg_mat,3)], 'Datatype', 'double', 'ChunkSize', [1, size(bg_mat,2), size(bg_mat,3)], 'Deflate', 9);
                h5write(H5filename, bg_name, bg_mat, [1 1 1], [size(bg_mat,1) size(bg_mat,2) size(bg_mat,3)]);
                
            end
            
            clear cellpath_mat bg_mat sisterList_mat
            
            if exist([maskOUTfilename '.mat'],'file')
                load(maskOUTfilename,'selected_cells');
                
                if exist(H5filename,'file')
                    fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                    if ~H5L.exists(fid,datasetname,'H5P_DEFAULT')
                        H5F.close(fid);
                        display(['Initializing ' H5filename ':' datasetname]);
                    else
                        H5L.delete(fid,datasetname,'H5P_DEFAULT');
                        display(['Overwriting ' H5filename ':' datasetname]);
                        H5F.close(fid);
                    end
                end
                
                h5create(H5filename, datasetname, [size(cellpath{last_tp},1), last_tp, 3, cellsize*2+1, cellsize*2+1], 'Datatype', 'uint8', 'ChunkSize', [1, 1, 3, cellsize*2+1, cellsize*2+1], 'Deflate', 9);
                
                for cellNo = selected_cells'
                    
                    newsegments = zeros(1,last_tp,3,cellsize*2+1,cellsize*2+1, 'uint8');
                    
                    nucName = ['nucmask_cell' num2str(cellNo)];
                    cytoName = ['cytomask_cell' num2str(cellNo)];
                    cellName = ['cellmask_cell' num2str(cellNo)];
                    
                    load(maskOUTfilename,nucName);
                    load(maskOUTfilename,cytoName);
                    load(maskOUTfilename,cellName);
                    
                    eval(['nucmask=' nucName ';']);
                    eval(['cytomask=' cytoName ';']);
                    eval(['cellmask=' cellName ';']);
                    
                    clear(nucName);
                    clear(cytoName);
                    clear(cellName);
                    
                    for tp=first_tp:last_tp
                        newsegments(1,tp,1,:,:) = nucmask(:,:,tp);
                        newsegments(1,tp,2,:,:) = cellmask(:,:,tp);
                        newsegments(1,tp,3,:,:) = cytomask(:,:,tp);
                    end
                    
                    clear nucmask;
                    clear cytomask;
                    clear cellmask;
                    
                    h5write(H5filename, datasetname, newsegments, [cellNo 1 1 1 1], [1 last_tp-first_tp+1 3 cellsize*2+1 cellsize*2+1]);
                end
                
                fid = H5F.open(H5filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                if ~H5L.exists(fid,selectedcells_name,'H5P_DEFAULT')
                    H5F.close(fid);
                    display(['Initializing ' H5filename ':' selectedcells_name]);
                else
                    H5L.delete(fid,selectedcells_name,'H5P_DEFAULT');
                    display(['Overwriting ' H5filename ':' selectedcells_name]);
                    H5F.close(fid);
                end
                
                h5create(H5filename, selectedcells_name, length(selected_cells), 'Datatype', 'uint32');
                h5write(H5filename, selectedcells_name, uint32(selected_cells));
                
            end
            
        end
    end
end