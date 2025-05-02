function out = f_gen(T,X,Q,N,M,L)
    % Notes: If X is passed without the rows for baseline and contrast,
    % X=X(1+2*Q:2*Q+M*L), then we don't need to pass M or L. If we replace 
    % reshape with squeeze, then we don't need Q or N. Reshape seems faster
    % than squeeze here, but reducing the number of passed parameters might
    % be better in total.
    out = reshape(pagemtimes(T,X(1+2*Q:2*Q+M*L)),Q,N);
end