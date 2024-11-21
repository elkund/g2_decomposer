function out = f_gen(T,X,Q,N,M,L)
    out = reshape(pagemtimes(T,X(1+2*Q:2*Q+M*L)),Q,N);
end