function result = KWW_compile_result(sample,series,X,q_value,t,q_powers,q_norm,t_norm,g2,g2_error)
% Turn a solution X into a human-readable "result" struct
% Q q bins, L components, M points per component

Q = length(q_value);
N = length(t);
L = length(q_powers);

baseline = X(1:Q);
contrast = X(1+Q:2*Q);

M=5E2; %set number of points in inverse transform
fit_eval_N = 100;
fit_eval_t = logspace(log10(t(1))-1,log10(t(end))+1,fit_eval_N); %extended delay vector

A_slice = 1+2*Q:2*Q+L-1;
G_slice = 1+2*Q+L-1:2*Q+L-1+L;

comp_strength = [1-sum(X(A_slice));X(A_slice)];

A = [1-sum(X(A_slice));X(A_slice)];
G = X(G_slice);
g = X(G_slice+L);
f_fits = f_gen_KWW(A,G,g,q_powers,q_value/q_norm,fit_eval_t/t_norm);
g2_fits = g2_gen(X,f_fits,Q);

f_fits_rough = f_gen_KWW(A,G,g,q_powers,q_value/q_norm,t/t_norm);
g2_fits_rough = g2_gen(X,f_fits_rough,Q);

f_comp_fits = zeros(Q,fit_eval_N,L);
for l = 1:L
    f_comp_fits(:,:,l) = f_gen_KWW(A(l),G(l),g(l),q_powers(l),q_value/q_norm,fit_eval_t/t_norm);
end

s = zeros(L,M);
component_dist = zeros(M,L);
t_powers = ceil(g)';

qn = q_value/q_norm;
tn = t/t_norm;

for l = 1:L
    tmp_poll = @(x) pollard(x,g(l)/t_powers(l) );
    mode = fminbnd( @(x) -tmp_poll(x) ,1E-2,1E4); %find mode
    max_val = tmp_poll(mode);

    % Define inverse transform limits by percentage of the maximum value.
    % (usually works fine)
    lo = fminbnd(@(x) (1E-3*max_val - tmp_poll(x))^2,1E-2,mode); %find low bound
    hi = fminbnd(@(x) (1E-3*max_val - tmp_poll(x))^2,mode,1E4); %find hi bound

    x = logspace(log10(lo),log10(hi),M);

    %change variable to rate
    s(l,:) = G(l)*x.^(1/t_powers(l));
    dxds = t_powers(l)*s.^(t_powers(l)-1)/(G(l)^t_powers(l));
    
    y = tmp_poll(x).*dxds;

    norm = trapz(s(l,:),y); % Normalization needed since we are missing a
    % little bit of the pdf.

    component_dist(:,l) = y/norm;
end

result= struct();

result.L = L;
result.Q = Q;
result.M = M;
result.N = N;
result.f_comp_fits = f_comp_fits;
result.baseline = baseline;
result.contrast = contrast;
result.G = G;
result.g = g;
result.component_strength = comp_strength;
result.component_dist = component_dist;
result.sample = sample;
result.series = series;
result.q_value =q_value;
result.s = s;
result.q_powers = q_powers*t_powers;
result.t_powers = t_powers';
result.g2 = g2;
result.g2_error = g2_error;
result.t = t;
result.fit_eval_t = fit_eval_t;
result.q_norm = q_norm;
result.t_norm = t_norm;
result.g2_fits = g2_fits;
result.g2_fits_rough = g2_fits_rough;
result.f_fits = f_fits;

end