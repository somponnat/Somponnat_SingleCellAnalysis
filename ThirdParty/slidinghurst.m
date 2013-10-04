function [output1, output2] = slidinghurst(signal,window)
for t=1:length(signal)-length(window)+1
    Hest = wfbmesti(signal(t+window-1));
    output1(t) = Hest(1);
    output2(t) = Hest(2);
end