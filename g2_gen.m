function out = g2_gen(X,fs,Q)
    %baselines X(1:Q), contrasts X(1+Q:2*Q), intermediate scattering
    %function f
    out = 1+X(1:Q)+X(1+Q:2*Q).*fs.^2;
end