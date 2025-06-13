function result = dist_fit_compile_result(sample,series,sol,q_value,t,q_powers,t_powers,q_norm,dist_strs,dist_Nx,dist_Nparams,lb,ub,t_norm,g2,g2_error)

Q = size(g2,1);
N = size(g2,2);
L = length(t_powers);

fit_Ndelays = 1E2; %number of delay times in g2,f reconstruction
fit_eval_t = logspace(log10(min(t))-1,log10(max(t))+1,fit_Ndelays);

dist_amplitudes = sol(1+2*Q:2*Q+L-1);
dist_amplitudes = [1-sum(dist_amplitudes) dist_amplitudes']';

dist_params = sol(1+2*Q+L-1:2*Q+L-1+sum(dist_Nparams));

f_fits = f_gen_from_dists(dist_amplitudes,dist_params,q_powers,t_powers,dist_strs,dist_Nx,dist_Nparams,fit_eval_t/t_norm,q_value/q_norm,Q,fit_Ndelays,L);

g2_fits = g2_gen(sol,f_fits,Q);

f_fits_rough = f_gen_from_dists(dist_amplitudes,dist_params,q_powers,t_powers,dist_strs,dist_Nx,dist_Nparams,t/t_norm,q_value/q_norm,Q,N,L);
g2_fits_rough = g2_gen(sol,f_fits_rough,Q);

f_comp_fits = zeros(Q,fit_Ndelays,L);

for l=1:L
    params_start_ind = 1+2*Q+L-1+sum(dist_Nparams(1:l-1));
    params_end_ind = 2*Q+L-1+sum(dist_Nparams(1:l));
    partmp = sol(params_start_ind:params_end_ind);

    %solsub = [sol(1:2*Q+L-1)' ']'
    f_comp_fits(:,:,l) = dist_amplitudes(l)*f_gen_from_dists(1,partmp,q_powers(l),t_powers(l),dist_strs(l),dist_Nx,dist_Nparams(l),fit_eval_t/t_norm,q_value/q_norm,Q,fit_Ndelays,1);
end

basefits = sol(1:Q);
contrastfits = sol(1+Q:2*Q);


% compile result as struct

result = struct();

result.sample = sample;
result.series = series;

M=dist_Nx; % number of points in pdf evaluation

comps = zeros(M,L);
s = zeros(1,M*L);


result.method = 'dist fit';
result.L = L;
result.Q = Q;
result.q_value = q_value;
result.M = M;
result.t = t;
result.g2 = g2;
result.g2_error = g2_error;
result.q_powers = q_powers;
result.t_powers = t_powers;
result.dist_strings = dist_strs;
result.dist_Nparams = dist_Nparams;
result.lo_bounds = lb;
result.hi_bounds = ub;
%result.init_guess = X0;
result.solution = sol;
result.dists = {};
result.component_strength = dist_amplitudes;
result.fit_eval_t = fit_eval_t;
result.contrast = contrastfits;
result.baseline = basefits;
result.g2_fits = g2_fits;
result.g2_fits_rough = g2_fits_rough;
result.f_fits = f_fits;
result.f_comp_fits = f_comp_fits;
result.q_norm = q_norm;
result.t_norm= t_norm;


for l = 1:L
    %find params for component l
    params_start_ind = 1+sum(dist_Nparams(1:l-1));
    params_end_ind = sum(dist_Nparams(1:l));

    makedist_args = num2cell(dist_params(params_start_ind:params_end_ind));
    dist = makedist(dist_strs{l},makedist_args{:});
    result.dists{l} = dist;

end

for l=1:L

    % this is specific to log normal dists, may not make sense in the general case
    % 
    % tot_np = sum(dist_Nparams); %total number of distribution parameters
    % 
    % current_ind = sum(dist_Nparams(1:l-1));
    % lo = exp(lb(end-tot_np+current_ind+1));   %plot between bounds on median
    % hi = exp(ub(end-tot_np+current_ind+1));   %plot between bounds on median

    %    params_start_ind = 1+sum(dist_Nparams(1:l-1));
    %params_end_ind = sum(dist_Nparams(1:l));

    %dist_args = num2cell(dist_params(params_start_ind:params_end_ind));
    dist = result.dists{l};%makedist(dist_strs{l},dist_args{:}); %These should be prepared and passed around instead of generated on each call.

    init=dist.median;

    l_prctl = fzero( @(x) cdf(dist,x)-1E-4,init); %Find the lower limit
    h_prctl = fzero( @(x) cdf(dist,x)-0.9999,init); %Find the upper limit
    x = logspace(log10(l_prctl),log10(h_prctl),dist_Nx);

    y = pdf(dist,x)';
    norm = trapz(x,y); % Normalization needed since we are missing a little bit of the pdf

    %x = logspace(log10(lo),log10(hi),M);

    comps(:,l)=y/norm;
    s(1+M*(l-1):M*l) = x;
end

result.component_dist = comps;
result.s = s;
end