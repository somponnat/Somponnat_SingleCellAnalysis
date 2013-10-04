function [peaktime,peakvalue,peakheight,Settlingheight,changedheight] = stimulationParams(t,y)

decompTS = awt1D(y);
currY = decompTS(:,4);

%STEP = stepinfo(currY,t,min(currY(length(currY)-10:end)));

%RiseTime = STEP.RiseTime;
%SettlingMin = STEP.SettlingMin;
%SettlingMax = STEP.SettlingMax;
%Overshoot = STEP.Overshoot;
peakvalue = max(currY);
peaktime = t(peakvalue==currY);
peakheight = peakvalue - min(currY(1:find(peakvalue==currY)));
Settlingheight = peakvalue - min(currY(find(peakvalue==currY):end));

currY = decompTS(:,5);
changedheight = currY(end) - currY(1);




% [outx,outy] =  slidingkendall(t,currY,1:15);
% 
% peakind  = find(min(abs(outy)) == abs(outy));
% peaktime = t(peakind);
% peakvalue = y(peakind);