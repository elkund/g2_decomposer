function [p,chisq, out_sol] = compute_GDF(g2,g2_error,X0s,opt_problem,Q,N,M,L,T,w,n_perts,D,lm)
% computes the loss in degrees of freedom according to Ye's generalized
% degrees of freedom. 

n_start = size(X0s,2);
n_jobs = n_start*(1+n_perts);

perts = normrnd(0,0.6*repmat(g2_error,1,1,1+n_perts),Q,N,1+n_perts);
perts(:,:,1) = 0;

par_iterations = [n_start,1+n_perts];
par_sols = zeros(2*Q+L*M,n_jobs); %temporary hold all solutions, including
% "bad" guesses for every perturbation
par_objvals = zeros(1,n_jobs);


parfor k = 1:n_jobs
    [guess,pert] = ind2sub(par_iterations,k);
    obj_fun  = @(x)obj_g2(T,x,g2+perts(:,:,pert),g2_error,Q,N,M,L,lm,w,D);
    par_prob = opt_problem;
    par_prob.objective = obj_fun;
    par_prob.X0 = X0s(:,guess);

    [tmpsol,tmpval,exitflag,~] = fmincon(par_prob);

    par_sols(:,k) = tmpsol;
    if exitflag > 0
        par_objvals(k) = tmpval;
    else
        par_objvals(k) = inf;
    end
    
end

pert_responses = zeros(Q,N,1+n_perts);
% collect result from best guess as perturbation responses
use_k =[];
for k=1:1+n_perts
    [val, ind] = min(par_objvals(1+(k-1)*n_start:k*n_start)); %Pick best run for each perturbation
    best_guess_sol = par_sols(:,(k-1)*n_start+ind);
 
    if val == inf
        if k ==1
            disp('Warning: No start points converged for unperturbed data. Skipping.')
        else
            disp('Warning: No start points converged for perturbation '+string(k-1)+', skipping.')
        end
    else
        pert_responses(:,:,k) = g2_gen(best_guess_sol,f_gen(T,best_guess_sol,Q,N,M,L),Q);
        use_k = [use_k k];
    end

    if k==1
        out_sol = best_guess_sol;
    end
end

%Compute loss in DOF, p
slopes = zeros(Q,N);
for q =1:Q
    for n=1:N
        x = squeeze(perts(q,n,use_k));
        y = squeeze(pert_responses(q,n,use_k));
        poly = polyfit(x,y,1);
        slopes(q,n) = poly(1);
    end
end
p = sum(slopes,'all');

chisq = RSS(pert_responses(:,:,1),g2,g2_error);

end