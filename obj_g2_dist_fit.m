function [out, grad] = obj_g2_dist_fit(X,q_deps,time_deps,dist_strs,dist_Nx,dist_Nparams,g2s,g2errs,qs,taus,Q,N,L)

% The following two lines do nothing if there is only one component. Can
% be better optimized?
dist_amplitudes = X(1+2*Q:2*Q+L-1);
dist_params = X(1+2*Q+L-1:2*Q+L-1+sum(dist_Nparams));

g2fits=g2_gen(X,f_gen_from_dists([1-sum(dist_amplitudes) dist_amplitudes']',dist_params,q_deps,time_deps,dist_strs,dist_Nx,dist_Nparams,taus,qs,Q,N,L),Q);

% Output to be minimized.
out = RSS(g2fits,g2s,g2errs);

end