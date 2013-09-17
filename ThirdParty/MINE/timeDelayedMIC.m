function [MIC,NLR] = timeDelayedMIC(data,varargin)

%This function calculates the MIC of a time series pair
%and its time delayed version
%
%Syntax : [MIC,NLR] = timeDelayedMIC(data,maxLag,style,exp,cv,c)
%
%INPUT:
%   data    : Input data(n of points,2)
%
%   maxLag  : maximum number of time delays (Lags);Default = 10
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
%
%       MIC - MIC(2*number of Lags + 1,1)
%       NLR - same size as MIC - Measure of nonlinearity
%
% See also: maxInformationCoef, autoMIC
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('data',@(x) isnumeric(x) && size(x,2) == 2);
ip.addParamValue('exp',0.5,@isscalar);
ip.addParamValue('cv',0.0,@isscalar);
ip.addParamValue('c',15,@isscalar);
nObs = size(data,1);
ip.addParamValue('maxLag',10,@(x) isscalar(x) && x <= round(nObs/4));


ip.parse(data,varargin{:});
maxLag = ip.Results.maxLag;
style  = '-allPairs';
exp    = ip.Results.exp;
cv     = ip.Results.cv;
c      = ip.Results.c;


MIC    = zeros(1,2*maxLag + 1);
NLR    = zeros(1,2*maxLag + 1);

[auxMIC,auxNLR] = maxInformationCoef([data(:,1) data(:,2)]',...
                        style,exp,cv,c);
MIC(1 + maxLag) = auxMIC(2,1);
NLR(1 + maxLag) = auxNLR(2,1);
    
for j = 1:maxLag
    %Shifting second variable into the past (shit that happens after the first variable)
    [auxMIC,auxNLR] = maxInformationCoef([data(1:end-j,1) data(j+1:end,2)]',...
                        style,exp,cv,c);
    
    MIC(j + maxLag) = auxMIC(2,1);
    NLR(j + maxLag) = auxNLR(2,1);
    
    %Shifting second variable into the future (shit that happens before the first variable)
    [auxMIC,auxNLR] = maxInformationCoef([data(j+1:end,1) data(1:end-j,2)]',...
                        style,exp,cv,c);
    
    MIC(j) = auxMIC(2,1);
    NLR(j) = auxNLR(2,1);
end
