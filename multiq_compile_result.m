function result = multiq_compile_result(sample,series,X,qs,delays,fit_eval_delays,s,T,fit_eval_T,w,q_deps,time_deps,g2s,g2errs,delay_norm,q_norm)
% Turn a solution X into a human-readable "result" struct
% Q q bins, L components, M points per component

Q = length(qs);
N = length(delays);
L = length(q_deps);
M = size(T,2)/L;

baseline_fits = X(1:Q);

contrast_fits = X(1+Q:2*Q);
components = reshape(X(1+2*Q:end),M,L);

comp_strengths = zeros(L,1);
for l=1:L
    comp_strengths(l) = w(1+M*(l-1):M*l)*components(:,l);
end

result= struct();

result.baseline_fits = baseline_fits;
result.contrast_fits = contrast_fits;
result.component_strengths = comp_strengths;
result.components = components;

result.sample = sample;
result.series = series;
result.qs =qs;
result.s = s;
result.T = T;
result.fit_eval_T = fit_eval_T;
result.q_deps = q_deps;
result.time_deps = time_deps;
result.g2s = g2s;
result.g2errs = g2errs;
result.delays = delays;
result.fit_eval_delays = fit_eval_delays;
result.q_norm = q_norm;
result.delay_norm = delay_norm;

fit_eval_N = length(result.fit_eval_delays);

f_eval_fits = f_gen(result.fit_eval_T,X,Q,fit_eval_N,M,L);
g2_eval_fits = g2_gen(X,f_eval_fits,Q);

result.g2_fits = g2_eval_fits;
result.f_fits = f_eval_fits;

f_comp_fits = zeros(Q,fit_eval_N,L);
for l = 1:L
    solsub = [X(1:2*Q)' X(1+2*Q+M*(l-1):2*Q+M*l)']';
    f_comp_fits(:,:,l) = f_gen(result.fit_eval_T(:,1+M*(l-1):M*l,:),solsub,Q,fit_eval_N,M,1);
end
result.f_comp_fits = f_comp_fits;

end