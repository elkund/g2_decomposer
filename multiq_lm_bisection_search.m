function [lm,sols,chisq,prob_reject] = multiq_lm_bisection_search(g2,g2_error,Q,N,M,L,T,w2,Delta_matrix,opt_problem,start_points,p_0,chisq_0,start_lm,max_its)
% Search for optimal regularizer weight with bisection method. The optimum
% is defined by Provencher's 50% rejection criterium.

thrsh_tolerance = 0.1;

sols = []; %holds the solutions for each iteration.
prob_reject = [];
lm = [];
chisq=[];
lm(1) = start_lm;
current_prob_reject = 0;

ms = MultiStart('UseParallel',true);

%Binary search limits, updated on each iteration
lm_a = 0;
lm_b = lm(1);

started = false;% set to true once lm such that prob_reject>0.5 is found
% (commence bisection search)

i=1;
while abs(current_prob_reject-0.5)>thrsh_tolerance && (i<=max_its)
    
    disp(' ')
    disp('Starting iteration ' + string(i) + ' with lm =' + string(lm(i)))

    obj_fun  = @(x)obj_g2(T,x,g2,g2_error,Q,N,M,L,lm(i),w2,Delta_matrix);
    opt_problem.objective = obj_fun;

    [Xsol,~,exitflag,~,~] = run(ms, opt_problem, start_points);

    if exitflag < 1
        disp('')
        disp('Warning: No start points converged.')
    end

    sols(:,i) = Xsol;
    g2_fits = g2_gen(Xsol,f_gen(T,Xsol,Q,N,M,L),Q);
    chisq(i) = RSS(g2_fits,g2,g2_error);

    current_prob_reject = fcdf((chisq(i)-chisq_0)/chisq_0,p_0,N*Q-p_0);
    prob_reject(i) = current_prob_reject;

    %start conditions
    if (~started) & (current_prob_reject<0.5)
        lm(i+1) = lm(i)*10; %Expand search area if initial run is smaller than 1
        lm_a= lm(i);
    elseif (~started) & (current_prob_reject>0.5+thrsh_tolerance)
        started = true;
    end

    % bisection method
    if started
        if current_prob_reject < 0.5
            lm_a = lm(i);
        else
            lm_b = lm(i);
        end
        lm(i+1) = (lm_a+lm_b)/2;
    end

    disp(' ')
    disp('Probability to reject is ' + string(current_prob_reject))

    i=i+1;
end

if i==max_its+1
    disp(' ')
    disp('Warning: probability  to reject = 0.5 was not reached in time')
end

disp(' ')
disp('Completed after ' +string(i-1) +' iterations with probability to reject '+ string(prob_reject(end)))

lm = lm(1:end-1);

end
