function [qs,delays,g2s,g2errs,ttc,ttc_q] = p10_load(year,prop,sample,series,q_sel,t_sel)
%Load g2s, and ttcs from P10 data
% year, prop and series are int, sample is a char array, e.g., 'sample_name'
% q_sel/t_sel are arrays of ints, e.g., [3 4 5]
% set q_sel/t_sel to 0 to select all
% Todo: use varargin and nargin to set default q_sel, t_sel

filename=strcat(sample,sprintf('_%05d_result.mat',series));
path = ['/asap3/petra3/gpfs/p10/' int2str(year) '/data/' int2str(prop) '/processed/RESULTS/' filename];
ccdimginfo = load(path).ccdimginfo;

all_qs = ccdimginfo.dmaskinfo(:,:,1)'; % Are these qs the the averages in each bin? Ask Michael.

g2s = double(squeeze(ccdimginfo.result.g2avg{1,1}));
g2errs = double(squeeze(ccdimginfo.result.g2avgErr{1,1})); %error bars
delays = ccdimginfo.result.delay{1,1};

if t_sel == 0
    t_sel = 1:length(delays);
end

if q_sel == 0
    q_sel = 1:length(all_qs);
end


qs = double(all_qs(q_sel));
g2s=g2s(q_sel,t_sel);
g2errs=g2errs(q_sel,t_sel);
delays=delays(t_sel);

ttc = ccdimginfo.result.twotime.C2t{1};
ttc_q = all_qs(ccdimginfo.result.twotime.bestq{1});

end