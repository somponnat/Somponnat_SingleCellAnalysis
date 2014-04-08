function [iMax, iWidth] = contrastValues(im,varargin)
% this function calculates mode and width of the intensity histogram

ip = inputParser;
ip.addRequired('im');
ip.addOptional('iStep', 1);
ip.parse(im, varargin{:});
iStep = ip.Results.iStep;

im = double(im);
imRange = min(im(:)):iStep:max(im(:));
nHist = hist(im(:), imRange);
nHist(1) = [];
imRange(1) = [];

nHist = filter([1 1 1], 1, nHist);
nMax = max(nHist);

% find the index of the peak
peakInd = find(nHist == max(nMax));
% find the index of the half peak from the lower end
halfPeakInd = max(find(nHist(1:peakInd) < 0.5 * nMax));
iMax = mean(imRange(peakInd));
iWidth = iMax - imRange(halfPeakInd);
end