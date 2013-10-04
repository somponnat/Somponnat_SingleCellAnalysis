function [outx outy] = slidingkendall(x,y,window)
for t=1:length(y)-length(window)+1
    outx(t) = x(t);
    outy(t) = corr(x(t+window-1),y(t+window-1),'type','Kendall');
    
end