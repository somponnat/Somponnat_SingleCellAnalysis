function [MIC,NLR,MAS,MEV,MCN] = maxInformationCoef(data,P,varargin)
% This is just a wrapper for the MINE java function 
%
%Synopsis: [MIC,NLR,MAS,MEV,MCN] = maxInformationCoef(data,style,exp,cv,c)
%
%INPUT:
%   data  : Input data(n of variables,n of points)
%
%   style : This option tells MINE which variable pairs to analyze. The value ?-allPairs? will cause MINE
%           to compare all pairs of variables against each other; ?-adjacentPairs? will compare consecutive pairs
%           of variables only; '-masterVariable i' will compare all variables only against the i-th variable;
%           ?-onePair i j? will compare only the i-th variable to the j-th variable; Variables are indexed from 0
%
%    cv : A floating point number indicating which percentage of the records need to have data in them for
%       both variables before those two variables are compared. Default value is 0.
%
%    exp : The exponent in the equation B(n) = n^exp . Default value is 0.6.
%
%    c : Determines by what factor clumps may outnumber columns when OptimizeXAxis is called. When
%       trying to partition the x-axis into x columns, the algorithm will start with at most cx clumps. Default
%       value is 15.
%
%OUTPUT:
%       All outputs have the same format (#Var,#Var)-  output between [var(i),var(j)] 
%       MIC - maximum information coeffient
%       NRL - Measure of nonlinearity (MIC - R^2)
%       MAS - Measure of monotonicity
%       MEV - Measure of complexity
%       MCN - Closeness of a functional form
%
% See also: autoMIC, timeDelayedMIC
%
%Reference: 
%Detecting Novel Associations in Large Data Sets
%David N. Reshef, ..., Pardis C. Sabeti
%Science 16 December 2011: 334 (6062), 1518-1524. [DOI:10.1126/science.1205438]
%
%Marco Vilela, 2012

%P = mfilename('fullpath');
%P(end-18:end) = [];

ip=inputParser;
ip.addRequired('data',@(x) isnumeric(x) && size(x,2)>2);
ip.addOptional('style','-allPairs',@ischar);
ip.addOptional('exp',0.6,@isscalar);
ip.addOptional('cv',0.0,@isscalar);
ip.addOptional('c',15,@isscalar);
ip.addOptional('norm',false,@isscalar);
ip.addOptional('prec',5,@(x)(isscalar(x) && isposint(x)));

ip.parse(data,varargin{:});
style = ip.Results.style;
exp   = ip.Results.exp;
cv    = ip.Results.cv;
c     = ip.Results.c;

%Normalization shouldn't effect correlation scores, and will reduce risk of
%rounding error when converting to .csv
if ip.Results.norm
    for j = 1:size(data,1)
        data(j,:) = data(j,:) - mean(data(j,:));
        data(j,:) = data(j,:) ./ std(data(j,:));
    end
end

[nVar,~] = size(data);
precN    = ip.Results.prec;
dlm      = ',';
MIC      = eye(nVar);%Maximum Information Coefficent
NLR      = nan(nVar);%Nonlinear measure
MAS      = nan(nVar);%Maximum Asymmetric Score
MEV      = nan(nVar);%Maximum Edge Value
MCN      = nan(nVar);%Minimum Cell Number

varName  = repmat('Var',nVar,1);
varName  = strcat(varName,num2str([1:nVar]'),repmat(dlm,nVar,1));
fileName = ['tmpFile_' num2str(feature('getpid')) '.csv'];
fid      = fopen(fileName, 'w');
format   = sprintf('%%.%dg%s',precN,dlm);

for i=1:nVar
    fprintf(fid, '%s',varName(i,:));
    str = sprintf(format,data(i,:));
    str = str(1:end-1);
    fwrite(fid, str, 'uchar');
    fprintf(fid, '\n');
end
fclose(fid);

comd = ['java -jar ',fullfile(P,'MINE.jar'),' ',fileName,' ',...
              style,...
              ' 0',...
              ' cv=',num2str(cv),...
              ' exp=',num2str(exp),...
              ' c=',num2str(c) ];
           
system(comd);
if ~cv
    sCv = '0.0';
        fStr = 'adjacentpairs';
else
    sCv = num2str(cv);
end

switch style(2:4)
    case('all')
        fStr = 'allpairs';
    case('adj')
        fStr = 'adjacentpairs';
    case('mas')
        fStr = ['mv=',style(end)];
    case('one')
        fStr = [ style(end-2),'-vs-',style(end) ];
end

        
    
rFile = strcat(fileName,',',fStr,',cv=',sCv,',B=n^',num2str(exp),...
               ',Results.csv');
           
javaResult = csv2cell(rFile);
% if ispc
%     system(['del ',' ',fileName(1:end-3),'*']);
% else
%     system(['rm ',' ',fileName(1:end-3),'*']);
% end

idx = cellfun(@(x) str2num(x(end)), javaResult(2:end,1:2) );

%_________________________________Output_______________________________________
MIC( sub2ind(size(MIC),idx(:,1),idx(:,2)) ) = cell2mat( javaResult(2:end,3) );
MIC( sub2ind(size(MIC),idx(:,2),idx(:,1)) ) = cell2mat( javaResult(2:end,3) );

%natural measure of nonlinearity by MIC ? r2
NLR( sub2ind(size(NLR),idx(:,1),idx(:,2)) ) = cell2mat( javaResult(2:end,4) );
NLR( sub2ind(size(NLR),idx(:,2),idx(:,1)) ) = cell2mat( javaResult(2:end,4) );

%detect deviation from monotonicity with the maximum asymmetry score
MAS( sub2ind(size(MAS),idx(:,1),idx(:,2)) ) = cell2mat( javaResult(2:end,5) );
MAS( sub2ind(size(MAS),idx(:,2),idx(:,1)) ) = cell2mat( javaResult(2:end,5) );

%complexity?
MEV( sub2ind(size(MEV),idx(:,1),idx(:,2)) ) = cell2mat( javaResult(2:end,6) );
MEV( sub2ind(size(MEV),idx(:,2),idx(:,1)) ) = cell2mat( javaResult(2:end,6) );

%and closeness to being a function
MCN( sub2ind(size(MCN),idx(:,1),idx(:,2)) ) = cell2mat( javaResult(2:end,7) );
MCN( sub2ind(size(MCN),idx(:,2),idx(:,1)) ) = cell2mat( javaResult(2:end,7) );
