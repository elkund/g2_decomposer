function f_out = f_gen_from_dists(dist_amplitudes,dist_params,q_deps,time_deps,dist_strs,dist_Nx,dist_Nparams,taus,qs,Q,N,L)

f_out = zeros(Q,N);

for l=1:L
    params_start_ind = 1+sum(dist_Nparams(1:l-1));
    params_end_ind = sum(dist_Nparams(1:l));

    dist_args = num2cell(dist_params(params_start_ind:params_end_ind));
    dist = makedist(dist_strs{l},dist_args{:}); %These should be prepared and passed around instead of generated on each call.

    init=dist.median;

    l_prctl = fzero( @(x) cdf(dist,x)-1E-4,init); %Find lower limit
    h_prctl = fzero( @(x) cdf(dist,x)-0.9999,init); %Find upper limit
    dist_x = logspace(log10(l_prctl),log10(h_prctl),dist_Nx);

    dist_y = pdf(dist,dist_x);
    norm = trapz(dist_x,dist_y); % Normalization needed since we are missing part of the pdf


    trapz_y = dist_y'.*exp(-(dist_x').^time_deps(l).*reshape(qs.^q_deps(l)'.*taus.^time_deps(l),1,Q,N));
    f_out = f_out + dist_amplitudes(l)*squeeze(trapz(dist_x,trapz_y))/norm;
end

end