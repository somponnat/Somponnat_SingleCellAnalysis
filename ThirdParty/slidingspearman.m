function [output] = slidingspearman(x,y,window)
for t=1:length(y)-length(window)+1
    output(t) = corr(x(t+window-1),y(t+window-1),'type','Kendall');

end