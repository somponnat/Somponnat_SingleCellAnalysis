function [aMIC,aNLR] = autoMIC(TS,varargin)

%This function calculates the maximum information coeficient on a given
% time series in relation to its time delayed version
%
%Synopsis : [aMIC,aNLR] = autoMIC(TS,style,exp,cv,c)
%
%INPUT:
%   TS    : Input data(n of points,n of variables)
%
%   maxLag  : maximum number of time delays (Lags); Default = 10
%
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
%       MIC - MIC(number of Lags,1)
%       NLR - same size as MIC - Measure of nonlinearity
%
% See also: maxInformationCoef, timeDelayedMIC
%
%Marco Vilela, 2012

ip=inputParser;
ip.addRequired('TS',@(x) isnumeric(x));
ip.addParamValue('exp',0.5,@isscalar);
ip.addParamValue('cv',0.0,@isscalar);
ip.addParamValue('c',15,@isscalar);
[nObs,nVar] = size(TS);
ip.addParamValue('maxLag',10,@(x) isscalar(x) && x <= round(nObs/4));

ip.parse(TS,varargin{:});
maxLag = ip.Results.maxLag;
style  = '-allPairs';
exp    = ip.Results.exp;
cv     = ip.Results.cv;
c      = ip.Results.c;

aMIC    = zeros(1,maxLag);
aNLR    = zeros(1,maxLag);
aMIC(1) = 1;

for i = 1:nVar
    for j = 2:maxLag
        [auxMIC,auxNLR] = maxInformationCoef([TS(1:end-j+1,i) TS(j:end,i)]',...
                            style,exp,cv,c);
                        
        aMIC(i,j) = auxMIC(2,1);
        aNLR(i,j) = auxNLR(2,1);
    end
end




