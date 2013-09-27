function [params,paramNames] = nonlinearTrendParams(x,y,varargin)
%NONLINEARTRENDPARAMS time series descriptors of nonlinear trends
%
% [params,paramNames] = nonlinearTrendParams(x,y)
% [params,paramNames] = nonlinearTrendParams(x,y,'OptionName',optionValue,...)
%
%
% Hunter Elliott
% 9/2013


%% ------------- Input --------------- %%

if nargin < 2 || isempty(x) || isempty(y) || numel(x) ~= numel(y)
    error('X and Y must be input as vectors of equal length!')
end

%Make sure we have column vectors
x = x(:);
y = y(:);

ip = inputParser;
ip.addParamValue('RunMIC',true,@islogical)
ip.addParamValue('DecomposeTimescales',true,@islogical)
ip.parse(varargin{:});

p = ip.Results;


%% ------------ Parameters ----------- %%
                                             
%String for naming time scales.              
tScaleBase = 'Timescale';

%Total number of parameters, for initializing arrays. 
if p.RunMIC
    nPar = 9;
else
    nPar = 4;
end

%% ----------- Decompostion --------- %%


if p.DecomposeTimescales
    %Do wavelet based timescale decomposition
    decompTS = awt1D(y);
    %And we include the raw timeseries as well, because some of the
    %descriptors will perform poorly on the decomposed timeseries
    decompTS = [y decompTS];
    
else
    decompTS = y;
    
end

%Create the strings for naming each timescale
nScales = size(decompTS,2);
scaleStr = {'Raw Data'};
if nScales > 1
    scaleStr = [scaleStr arrayfun(@(x)([tScaleBase ' ' num2str(x)]),1:nScales-1,'Unif',0)];
end


%% ------------ Processing ---------- %%

nTot = nScales*nPar;
paramNames = cell(nTot,1);
params = nan(nTot,1);

for iScale = 1:nScales
           
    currY = decompTS(:,iScale);
    
    %We just hard-code the function numbers and order as using function
    %handle arrays adds unnecessary complexity and / or redundant function
    %calls            
    
    Hest = wfbmesti(currY);
    
    iPar = (iScale-1)*nPar + 1;
    paramNames{iPar} = [scaleStr{iScale} ' Hurst Exponent Estimate 1'];
    params(iPar) = Hest(1);
    
    iPar = (iScale-1)*nPar + 2;
    paramNames{iPar} = [scaleStr{iScale} ' Hurst Exponent Estimate 2'];
    params(iPar) = Hest(2);
            
    
    
    iPar = (iScale-1)*nPar + 3;
    paramNames{iPar} = [scaleStr{iScale} ' Kendall''s Tau'];
    params(iPar) = corr(x,currY,'type','Kendall');
    
    iPar = (iScale-1)*nPar + 4;
    paramNames{iPar} = [scaleStr{iScale} ' Spearman''s Rho'];
    params(iPar) = corr(x,currY,'type','Spearman');
    
    if p.RunMIC
    

        [MIC,NLR,MAS,MEV,MCN] = maxInformationCoef(vertcat(x',currY'),[pwd filesep 'ThirdParty' filesep 'MINE'],'-onePair 0 1','prec',12,'norm',true);
        %[MIC,NLR,MAS,MEV,MCN] = maxInformationCoef(vertcat(x',currY'),'-onePair 0 1');

        iPar = (iScale-1)*nPar + 5;
        paramNames{iPar} = [scaleStr{iScale} ' Maximum Information Coefficient'];
        params(iPar) = MIC(2,1);

        iPar = (iScale-1)*nPar + 6;
        paramNames{iPar} = [scaleStr{iScale} ' Nonlinearity'];
        params(iPar) = NLR(2,1);

        iPar = (iScale-1)*nPar + 7;
        paramNames{iPar} = [scaleStr{iScale} ' Monotonicity'];
        params(iPar) = MAS(2,1);

        iPar = (iScale-1)*nPar + 8;
        paramNames{iPar} = [scaleStr{iScale} ' Complexity'];
        params(iPar) = MEV(2,1);

        iPar = (iScale-1)*nPar + 9;
        paramNames{iPar} = [scaleStr{iScale} ' Functionality'];
        params(iPar) = MCN(2,1);
    end    
   
end




