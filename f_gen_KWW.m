function out = f_gen_KWW(A,G,g,q_power_law,q_value,t)
% Multi-component, multi-q, multi-t KWW function.
out = squeeze(pagemtimes(exp(-(G'.*q_value'.^(q_power_law).*reshape(t,1,1,[])).^(g')),A));
end