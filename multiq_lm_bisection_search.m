function [result] = multiq_lm_bisection_search(g2s,g2errs,X0s,T,s,w,lb,ub,Aeq,beq,opts,start_lm,max_its)
% Search for optimal lagrange multiplier  with bisection method
% Optimal is defined as reduced chi sq ~ 1
% Degrees of freedom computed with Ye's method for
% generalized degrees of freedom
% Uses parfor, parallell pool should be set up before call to this function

% This script is a big mess. It works, but needs to be rewritten.

Q=size(g2s,1);
N=size(g2s,2);
L=size(s,1);
M=size(s,2);

n_start = size(X0s,2);
n_perts = round(Q*N/7);
n_jobs = n_start*(1+n_perts);

%max_its=1;
redchi_thrs = 0.05;

redchisq = 0;

%all_sols = zeros(X,n_jobs,max_its);
sols = zeros(2*Q+L*M,max_its+1); %holds the unperturbed solutions for each iteration. Leave first page as zero
redchisqs = zeros(1,max_its+1); %set redchisq for lm=0 to zero, we don't actually have to compute it as it's sure to be smaller than 1.
DOFs = zeros(1,max_its+1);
lms = zeros(1,max_its+2); %Keep "first" lm 0, one extra column to hold last lm
lms(:,2) = start_lm;%9E-4;%/2^17; %we will start at iteration 2

%Binary search limits, updated on each iteration
lm_a = 0;%zeros(Q,1);
lm_b = lms(:,2);%zeros(Q,1);
sign_a = -1;%ones(Q,1);

termination_i = max_its+1;
started = false;% set to true once lm such that redchisq>1 is found (commence bisection)

pert_responses = zeros(Q,N,1+n_perts);
best_guesses = zeros(1,1+max_its);

par_iterations = [n_start,1+n_perts];
par_sols = zeros(2*Q+L*M,n_jobs); %temporary hold all solutions, "bad"
%  guesses and perturbation
par_objvals = zeros(1,n_jobs);


i=2; % start at two, no need to carry out the fit at lm =0
while abs(redchisq-1)>redchi_thrs & (i<max_its+2)

    lm = lms(i);
    'Starting iteration ' + string(i) + ' with lm =' + string(lm)

    perts = normrnd(0,0.6*repmat(g2errs,1,1,1+n_perts),Q,N,1+n_perts); %generate new perts for each iteration
    perts(:,:,1) = 0;

    parfor k = 1:n_jobs
        [guess,pert] = ind2sub(par_iterations,k);
        obj_fun  = @(x)obj_g2(T,w,x,g2s+perts(:,:,pert),g2errs,Q,N,s,M,L,lm);
        
        [tmpsol,tmpval] = fmincon(createOptimProblem('fmincon','objective',obj_fun,'x0',X0s(:,guess),'lb',lb,'ub',ub,'Aeq',Aeq,'beq',beq,'options',opts));

        par_sols(:,k) = tmpsol;
        par_objvals(k) = tmpval;
        %all_sols(:,k,i) = tmpsol;
        
    end

    % collect result from best guess as perturbation responses
    for k=1:1+n_perts
        [~, ind] = min(par_objvals(1+(k-1)*n_start:k*n_start));
        best_guesses(i) = ind;
        best_guess_sol = par_sols(:,(k-1)*n_start+ind);
        pert_responses(:,:,k) = g2_gen(best_guess_sol,g1_gen(T,best_guess_sol,Q,N,M,L),Q);
        if k==1
            sols(:,i) = best_guess_sol;
        end
    end

    %Compute red chi sq
    slopes = zeros(Q,N);
    for q =1:Q
        for n=1:N
            x = squeeze(perts(q,n,:));
            y = squeeze(pert_responses(q,n,:));
            poly = polyfit(x,y,1);
            slopes(q,n) = poly(1);
        end
    end
    
    nu = sum(slopes,'all');
    DOF = Q*N - nu;

    msd = MSD(pert_responses(:,:,1),g2s,g2errs); %MSD ~ DOF
    redchisq = msd/DOF; %~1

    % bisection method
    if (~started) & (redchisq<1-redchi_thrs)
        lms(i+1) = lm*10; %Expand search area if initial run is smaller than 1
        lm_a= lm;
    elseif (~started) & (redchisq>1+redchi_thrs)
        started = true;
    elseif abs(redchisq-1)<redchi_thrs
        termination_i = i;
        'Completed after ' +string(i) +' iterations with reduced chi square '+ string(redchisq)
    end

    if started 
        if sign(redchisq-1) == sign_a
            lm_a = lm;
        else
            lm_b = lm;
        end
        lms(i+1) = (lm_a+lm_b)/2;
    end

    'Reduced chi squared is ' + string(redchisq)
    'DOF is ' + string(DOF)
    
    %assign values to keep
    DOFs(i) = DOF;
    redchisqs(i) = redchisq;

    i=i+1;
end

if i==max_its+2
    'Warning, redchisq=1 was not reached in time'
end

%Xsol = sols(:,termination_i);

%base_fits = Xsol(1:Q);

%normalization  = w*Xsol(1+2*Q:end); %unnecesaary step since when norm=1 is required in the fitting
%contrast_fits = Xsol(1+Q:2*Q);%*normalization^2;
%comps = reshape(Xsol(1+2*Q:end),M,L);%/normalization;

%comp_strengths = zeros(L,1);
%for l=1:L
%    comp_strengths(l) = w(1+M*(l-1):M*l)*comps(:,l);
%    %comps(:,l) = comp_strengths(l)*comps(:,l);
%end

result= struct();

result.lms = lms(1:termination_i);
result.redchisqs = redchisqs(1:termination_i);
result.DOFs = DOFs(1:termination_i);
result.solutions = sols(:,1:termination_i);

%result.base_fits = base_fits;
%result.contrast_fits = contrast_fits;
%result.comp_strengths = comp_strengths;
%result.comps = comps;

end
