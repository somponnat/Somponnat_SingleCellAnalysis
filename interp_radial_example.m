close all
nsample = 20;

t = linspace(eps,2*pi,30);
% t = t(1:end-1);

xin = (1+.1*sin(6*t)).*cos(t)+1;
yin = (1+.1*sin(6*t)).*sin(t)+1;

x1in = .8*cos(t)+.5;
y1in = .8*sin(t);

plot(xin,yin,'o-',x1in,y1in,'o-');
hold on
% plot(xfine,yfine,x1fine,y1fine)

[interpx interpy] = interp_radial(xin,yin,x1in,y1in,nsample);
plot(interpx,interpy,'ro-')
