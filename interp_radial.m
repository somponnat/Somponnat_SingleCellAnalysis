function [interpx interpy] = interp_radial(xin,yin,x1in,y1in,nsample)
    offset_xin = mean(xin);
    offset_yin = mean(yin);
    offset_x1in = mean(x1in);
    offset_y1in = mean(y1in);

    xin = xin - offset_xin;
    yin = yin - offset_yin;
    x1in = x1in - offset_x1in;
    y1in = y1in - offset_y1in;

    be = atan(yin./xin);
    be(be<0) = pi+be(be<0);
    be(yin<0) = pi+be(yin<0);
    be1 = atan(y1in./x1in);
    be1(be1<0) = pi+be1(be1<0);
    be1(y1in<0) = pi+be1(y1in<0);

    tfine = linspace(0,2*pi,501);
    xfine = interp1(be,xin,tfine);
    yfine = interp1(be,yin,tfine);
    x1fine = interp1(be1,x1in,tfine);
    y1fine = interp1(be1,y1in,tfine);

    as = linspace(0,2*pi,nsample);

    interpx = nan(size(as));
    interpy = nan(size(as));
    for ai = 1:length(as)-1
        ind = tfine > as(ai) & tfine < as(ai+1);
        xm = mean(xfine(ind));
        ym = mean(yfine(ind));
        x1m = mean(x1fine(ind));
        y1m = mean(y1fine(ind));

        interpx(ai) = mean([xm x1m]);
        interpy(ai) = mean([ym y1m]);
    end

    interpx(end) = interpx(1);
    interpy(end) = interpy(1);

    interpx = interpx + mean([offset_xin offset_x1in]);
    interpy = interpy + mean([offset_yin offset_y1in]);
end