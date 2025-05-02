function out = g2_gen(X,f,Q)
    %baselines X(1:Q), contrasts X(1+Q:2*Q), intermediate scattering
    %function f. If we pass X(1:Q), and X(1+Q:2*Q) separately, we can skip
    % the Q argument and the indexing operation.
    out = 1+X(1:Q)+X(1+Q:2*Q).*f.^2;
end