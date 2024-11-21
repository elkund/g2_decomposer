function out = MSD(g2fits,g2s,g2errs)
    devs = (g2fits-g2s);
    out = sum((devs./g2errs).^2,[1 2]); %Why is this faster than sumsqr()??
end