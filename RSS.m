function out = RSS(g2fits,g2s,g2errs)
    devs = (g2s-g2fits);
    out = sum((devs./g2errs).^2,[1 2]); %Why is this faster than sumsqr()??
end