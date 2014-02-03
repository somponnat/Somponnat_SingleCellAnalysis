field = 1;
SourceF = 'Q:\sorger\data\NIC\Pat\07-04-2013';
for row=2:8
    for col = 3:10
        H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
        if exist(fullfile(SourceF,H5filename),'file')
            cleanH5file(row,col,field,SourceF)
        end
    end
end