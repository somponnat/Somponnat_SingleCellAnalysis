function parallelmoviegen()
rows = 1:8;
cols = 4:5;
fields = 1;
plane = 1;
templateCH = 1;
nominCH = -1;
denominCH = -1;
tps = [1 55];
lind = 1;

matlabpool open;

for row = rows
    for col = cols
        for field = fields
                loopparam(lind,:) = [row col field];
                lind = lind+1;
        end
    end
end


parfor lind = 1:size(loopparam,1)
    row = loopparam(lind,1);
    col = loopparam(lind,2);
    field = loopparam(lind,3);
    GenMov_commandline(row,col,field,plane,templateCH,nominCH,denominCH,tps)
end
matlabpool close;