function cleanlist_commandline(celltrackOUT)
% hObject    handle to pushbutton_cleanList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load(celltrackOUT);
first_tp = 1;
last_tp = length(cellpath);

if ~isempty(cellpath)
    

    sis_cellpath = cell(last_tp,1);
    sis_sisterList = cell(last_tp,1);
    % Determine cells with sisters
    withSisInd = find(sisterList{last_tp}(:,1)~=-1);
    
    if ~isempty(withSisInd)
        
        cInd=1;
        loopInd=1;
        lastInd(loopInd)=1;
        while ~isempty(withSisInd)
            firstSis = withSisInd(1);
            secondSis = sisterList{last_tp}(withSisInd(1),1);
            sis_gInd = find(sisterList{last_tp}(:,1)==firstSis | sisterList{last_tp}(:,1)==secondSis);
            
            for t = last_tp:-1:first_tp
                if  t~=last_tp
                    sis_sisterList{t} = -1*ones(size(sis_sisterList{last_tp}));
                else
                    
                    cInd=lastInd(loopInd);
                    
                    for s=1:length(sis_gInd)
                        sis_table{sis_gInd(s)} = s+cInd-1;
                    end
                    
                    for s=1:length(sis_gInd)
                        sis_cellpath{t}(cInd,:)   = cellpath{t}(sis_gInd(s),:);
                        oldList = sisterList{t}(sis_gInd(s),:);
                        posL = find(oldList ~= -1);
                        
                        if ~isempty(posL)
                            switch length(posL)
                                case 1
                                    sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} -1 -1];
                                case 2
                                    sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} sis_table{oldList(2)} -1];
                                case 3
                                    sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} sis_table{oldList(2)} sis_table{oldList(3)}];
                            end
                        end
                        cInd = cInd+1;
                    end
                    
                end
                
            end
            loopInd=loopInd+1;
            lastInd(loopInd) = cInd;
            withSisInd = setdiff(withSisInd,sis_gInd);
            clear sis_table;
            
        end
        
        cInd=1;
        loopInd=1;
        lastInd(loopInd)=1;
        withSisInd = find(sisterList{last_tp}(:,1)~=-1);
        while ~isempty(withSisInd)
            firstSis = withSisInd(1);
            secondSis = sisterList{last_tp}(withSisInd(1),1);
            sis_gInd = find(sisterList{last_tp}(:,1)==firstSis | sisterList{last_tp}(:,1)==secondSis);
            
            for t = last_tp:-1:first_tp
                
                cInd=lastInd(loopInd);
                
                for s=1:length(sis_gInd)
                    sis_table{sis_gInd(s)} = s+cInd-1;
                end
                
                for s=1:length(sis_gInd)
                    sis_cellpath{t}(cInd,:)   = cellpath{t}(sis_gInd(s),:);
                    oldList = sisterList{t}(sis_gInd(s),:);
                    posL = find(oldList ~= -1);
                    if ~isempty(posL)
                        switch length(posL)
                            case 1
                                sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} -1 -1];
                            case 2
                                sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} sis_table{oldList(2)} -1];
                            case 3
                                sis_sisterList{t}(cInd,:) = [sis_table{oldList(1)} sis_table{oldList(2)} sis_table{oldList(3)}];
                        end
                    end
                    cInd = cInd+1;
                end
                
            end
            loopInd=loopInd+1;
            lastInd(loopInd) = cInd;
            withSisInd = setdiff(withSisInd,sis_gInd);
            clear sis_table;
        end
        
    end
    
    % Determine cells without sisters
    noSisInd = find(sisterList{last_tp}(:,1)==-1 & cellpath{last_tp}(:,1)~=-1);
    for t = first_tp:last_tp
        nosis_cellpath{t}   = cellpath{t}(noSisInd,:);
        nosis_sisterList{t} = sisterList{t}(noSisInd,:) ;
    end
    
    [single_cellpath single_sisterList couple_cellpath couple_sisterList] = massWedding(nosis_cellpath,nosis_sisterList,first_tp,last_tp);

    if ~isempty(couple_cellpath)
        for t = first_tp:last_tp
            if ~isempty(sis_cellpath)
                oldsissize = size(sis_cellpath{t},1);
                new_cellpath{t} =   sis_cellpath{t};
                new_sisterList{t} = sis_sisterList{t};
            else
                oldsissize = 0;
            end
            
            addsissize = size(couple_cellpath{t},1);
            for c = 1 : addsissize
                new_cellpath{t}(oldsissize+c,:)   =  couple_cellpath{t}(c,:);
                old_sisters = couple_sisterList{t}(c,:);
                PosSis = find(old_sisters>0);
                new_sisterList{t}(oldsissize+c,:) =  old_sisters;
                if ~isempty(PosSis)
                    for ss=1:length(PosSis)
                        new_sisterList{t}(oldsissize+c,PosSis(ss)) = old_sisters(PosSis(ss))+oldsissize;
                    end
                end
            end
        end
        sis_cellpath = new_cellpath;
        sis_sisterList = new_sisterList;

    end
    
    for t = first_tp:last_tp
        new_cellpath{t} =[sis_cellpath{t};single_cellpath{t}];
        new_sisterList{t} = [sis_sisterList{t};single_sisterList{t}];
    end
    
    cellpath   = new_cellpath;
    sisterList = new_sisterList;
    assignin('base','cellpath',cellpath);
    save(celltrackOUT,'cellpath','sisterList','bg','-v7.3');
end


function [single_cellpath single_sisterList couple_cellpath couple_sisterList] = massWedding(nosis_cellpath,nosis_sisterList,first_tp,last_tp)
single_cellpath   = nosis_cellpath;
single_sisterList = nosis_sisterList;

couple_cellpath   = cell(last_tp,1);
couple_sisterList = cell(last_tp,1);

Rad = 3;
%create distance matrix for first time point
D = pdist(nosis_cellpath{first_tp},'euclidean');
disMat = squareform(D);
assignin('base','disMat',disMat);

noSisList = [];
SisList = [];
out_Ind = 1;
leftInd = 1:size(disMat,1)';
while ~isempty(leftInd)
    sis1Ind = find(disMat(leftInd(1),:)<Rad);
    sisInd = find(sis1Ind~=leftInd(1));
    if isempty(sisInd)
        noSisList = [noSisList;leftInd(1)];
        cList = [leftInd(1)];
    else
        switch length(sisInd)
            case 1
                ind_dist = zeros(1,last_tp);
                for t=first_tp:last_tp
                    ind_dist(t) = pdist([nosis_cellpath{t}(leftInd(1),:);nosis_cellpath{t}(sis1Ind(sisInd(1)),:)]);
                end
                diff_dist = diff(ind_dist);
                jumpInd = find(diff_dist>0,1,'first');
                if isempty(jumpInd)
                    splitF = 1;
                else
                    splitF = jumpInd+1;
                end
                for t=first_tp:splitF-1
                    couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                    couple_sisterList{t}(out_Ind+1,:) = [-1 -1 -1];
                    couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(leftInd(1),:);
                    couple_cellpath{t}(out_Ind+1,:) = [-1 -1];
                end 
                for t=splitF:last_tp
                    couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                    couple_sisterList{t}(out_Ind+1,:) = [out_Ind -1 -1];
                    couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(leftInd(1),:);
                    couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(sis1Ind(sisInd(1)),:);
                end
                SisList = [SisList;leftInd(1);sis1Ind(sisInd(1))];
                cList = [leftInd(1);sis1Ind(sisInd(1))];
                out_Ind=out_Ind+2;
            case 2
                ind_dist = zeros(1,last_tp);
                for t=first_tp:last_tp
                    ind_dist(t) = pdist([nosis_cellpath{t}(leftInd(1),:);nosis_cellpath{t}(sis1Ind(sisInd(1)),:)]);
                end
                diff_dist = diff(ind_dist);
                jumpInd = find(diff_dist>0,1,'first');
                if isempty(jumpInd)
                    splitF = 1;
                else
                    splitF = jumpInd+1;
                end
                
                for t=first_tp:splitF-1
                    couple_sisterList{t}(out_Ind,:) = [-1 -1 -1];
                    couple_sisterList{t}(out_Ind+1,:) = [-1 -1 -1];
                    couple_sisterList{t}(out_Ind+2,:) = [-1 -1 -1];
                    couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(leftInd(1),:);
                    couple_cellpath{t}(out_Ind+1,:) = [-1 -1];
                    couple_cellpath{t}(out_Ind+2,:) = [-1 -1];
                end
                for t=splitF:last_tp
                    couple_sisterList{t}(out_Ind,:) = [out_Ind+1 -1 -1];
                    couple_sisterList{t}(out_Ind+1,:) = [out_Ind -1 -1];
                    couple_sisterList{t}(out_Ind+2,:) = [out_Ind -1 -1];
                    couple_cellpath{t}(out_Ind,:) = nosis_cellpath{t}(leftInd(1),:);
                    couple_cellpath{t}(out_Ind+1,:) = nosis_cellpath{t}(sis1Ind(sisInd(1)),:);
                    couple_cellpath{t}(out_Ind+2,:) = nosis_cellpath{t}(sis1Ind(sisInd(2)),:);
                end
                SisList = [SisList;leftInd(1);sis1Ind(sisInd(1));sis1Ind(sisInd(2))];
                cList = [leftInd(1);sis1Ind(sisInd(1));sis1Ind(sisInd(2))];
                out_Ind=out_Ind+3;
            otherwise
                noSisList = [noSisList;leftInd(1)];
                cList = [leftInd(1)];
        end
    end
    leftInd = setdiff(leftInd,cList);
  
end

for t=first_tp:last_tp
    single_cellpath{t} = nosis_cellpath{t}(noSisList,:);
    single_sisterList{t} = nosis_sisterList{t}(noSisList,:);
end