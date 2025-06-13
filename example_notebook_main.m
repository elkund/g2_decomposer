%% Load example data and optimization settings

% Uncomment to load data and settings. Use the settings as template for new
% datasets.

%%%%%%%%%
% Colloidal particles at 313 K
%%%%%%%%
load("example_data/PRE2021_313K.mat")
% load("example_settings/PRE2021_313K_2comp.mat")
load("example_settings/PRE2021_313K_1comp.mat")


%%%%%%%%%
% Colloidal partices at 288 K
% %%%%%%%%
% load("example_data/PRE2021_288K.mat");
% load("example_settings/PRE2021_288K_2comp.mat");

%%%%%%%%%
% Amorphous ice at 125 K
%%%%%%%%%
% load("example_data/ASW_125K.mat")
% load("example_settings/ASW_125K.mat");

% Set array sizes
N=length(t);
Q=length(q_value);

%% Display data

% Settings for plotting g2s
set(groot,'defaulttextinterpreter','latex')
xlbl = '$\tau \quad [s]$';
ylbl = '$g_2$';

figure();

ax1=gca;
c_map=parula(Q+1);

hold on
for q = 1:Q
    label = strjoin([compose("%5.2e",q_value(q)') ' ' char(197) '^{-1}'],'');
    errorbar(t,g2(q,:),g2_error(q,:),'o', 'color',c_map(q,:),'DisplayName',label);
end

lgd=legend('Location','southwest');
ax1.XScale="log";
xlabel(xlbl);
ylabel(ylbl);


%% MULTIQ analysis setup

% Set dynamic component power laws. For example,
%   f = integral w. respect to s of
%       rho_1 e^(-s q^q_powers(1)*t^t_powers(1))
%      +rho_2 e^(-s q^q_powers(2)*t^t_powers(2))
%      + ...
% Program will solve for rho_1, rho_2, and so on.
% e.g.
% t_powers = [1 2];
% q_powers = [2 2];
% !! Already assigned during load!!

M = 150; %number of points in transformation, per component (total M*L)
L = length(q_powers);   %number of components

% Set the integration limits, per component. Shape of s_range should be L
% rows by 2 columns:
%      s_range=   [ lo 1, hi 1; 
%                   lo 2, hi 2;
%                   lo 3, hi 3;
%                   ... ]
% s_range = % already set during load.

% q and delay normalization. These numbers will be used to normalize the
% design matrix, e^(-(q/q_norm)^q_power*(s*t/t_norm)^t_power).
q_norm = max(q_value);
t_norm = max(t);

% Set up model arrays.
% Weighted kernel values T, quadrature weights w, integrand evaluation
% points s, differentitation matrix D (for regularizer)
[T,w,s,D] = multiq_setup_arrays(q_value,t,s_range,q_powers,t_powers,M,'q_norm',q_norm,'t_norm',t_norm);

% Set up another T array with more delay values. This is used to
% display extrapolated fit curves, and can be higher
% res (if fit_eval_N is larger than N). We could also interpolate/extrapolate along
% q by giving more q values.

% Number of delay points in back-transformed results
fit_eval_N=100;

fit_eval_t = logspace(log10(t(1))-1,log10(t(end))+1,fit_eval_N); %extended delay vector
[fit_eval_T,~,~,~] = multiq_setup_arrays(q_value,fit_eval_t,s_range,q_powers,t_powers,M,'q_norm',q_norm,'t_norm',t_norm);

%% Bounds and constraints on solution

% Equality constraints.
% Fmincon will find a solution X such that Aeq*X=beq.

reg_order=2; %Regularizer derivative order (normally 2)
Aeq = zeros(2,2*Q+M*L);

% Boundary condiitions (sum of values on boundary). We set the first and
% last n=reg_order to zero.
for l=1:L
    Aeq(2,1+2*Q+M*(l-1):1+2*Q+M*(l-1)+reg_order-1) = 1;
    Aeq(2,2*Q+M*l-reg_order+1:2*Q+M*l) = 1;
end

%Constrain integral of |f| to 1 
Aeq(1,1+2*Q:end) = w;

beq = [1 0]'; %total = 1, boundaries = 0

% Upper and lower bounds.
% fmincon will find a solution such that lb =< X =< ub

% Bounds on component values
lb = zeros(2*Q+M*L,1); %Minimum is zero.
ub = ones(2*Q+M*L,1);
ub(2*Q+1:end) = 1./w; %Every element contributes maximum area 1 to transformed function.

% Bounds on baselines. We use the same number for every q, but per q values
% can be supplied as well.
lb(1:Q) = baseline_lb; %baseline lower bound
ub(1:Q) = baseline_ub; %baseline uper bound

% bounds on constrast.
lb(Q+1:2*Q) = contrast_lb; %contrast lower bound 
ub(Q+1:2*Q) = contrast_ub; %contrast upper bound

%% Initial guesses

% Generate a list of initial guesees. Use n_start bump functions per
% component, total (n_start+1)^L-1 starting points. Baseline and contrast
% start at midpoint between bounds.
X0s = multiq_multistart_guess_gen(Q,M,n_start,L,w,lb,ub);

%% Initialize parallel pool

% Local parallel pool.
parpool(40);

% % On cluster
% c = parcluster;
% 
% %%%%%
% % On Maxwell, set up profile before. Instructions at:
% % https://confluence.desy.de/display/MXW/Getting-Started-with-Parallel-Computing-using-MATLAB-on-Maxwell_272733172.html
% %%%%%
% 
% % Specific to DESY's Maxwell cluster
% c.AdditionalProperties.QueueName = 'exfel';
% c.AdditionalProperties.WallTime = '08:00:00';
% c.saveProfile
% 
% c.parpool(40);

%However, we find that parpool(n) with "processes" profile usually performs
%better.

%% Set up multistart solve with fixed regularizer weight

% Objective function. Change lm to increase/decrease regularizer weight.
lm = start_lm;
%Starting at small lm gives an overfitted result. We will search for
%optinal lm later. You can try different values to see the effect on the
%result.

% Objective function
obj_fun  = @(x)obj_g2(T,x,g2,g2_error,Q,N,M,L,lm,w,D);

% Create problem formulation

% Optimization options
options = optimoptions(@fmincon,'Display','off','SpecifyObjectiveGradient',false,'ScaleProblem',true,'Algorithm','sqp','MaxFunctionEvaluations',length(X0s(:,1))*1E3,'MaxIterations',max_opt_iterations);
% Sqp alghorithm works well. Specify objective gradient seems to make most
% things worse.
% Scale problem helps, even though we normalized q and t.

prob = createOptimProblem('fmincon','objective',obj_fun,'x0',X0s(:,1),'lb',lb,'ub',ub,'Aeq',Aeq,'beq',beq,'options',options);

% Create "multistart" object
ms = MultiStart('UseParallel',true);
custom_start_points = CustomStartPointSet(X0s');

%% Run multistart
% find the best solution given the start points in X0s

tic

[Xsol1,fval,exitflag,output,all_solutions] = run(ms, prob, custom_start_points);

toc
%% Show result
transform_s = 1; %set to 0 to plot solution against dimless param s.
% Set to 1 to convert to input units.

% compile solution as human-readable struct
result_fixed_lm = multiq_compile_result(sample,series,Xsol1,q_value,t,fit_eval_t,s,T,fit_eval_T,w,q_powers,t_powers,g2,g2_error,t_norm,q_norm);

%plot
multiq_plot_result(result_fixed_lm,transform_s);


%% Compute DOF at small λ

%%%%%
% This calculation can take a long time. It is advised to start as many
% parallell workers as possible. Especially if solving for L>1. Decreasing
% n_perts or n_start shortens computation time. Hihger n_perts gives
% better precision for DOF, higher n_start increases chance of finding a
% 'close-to-global' optimum. Keep n_start consistent during analysis!
%
% The choice of reference λ can have an influence on the result for
% "best" λ (next section). Like Provencher, we have found that a small non-zero
% λ gives more reasonably consistent results.
%
% If you encounter warnings about solvers not converging, increase
% max_opt_iterations or find better initial guesses.
%%%%%

n_perts = 50;%fix(N*Q/4); %This number is a compromise. n_perts>N*Q is ideal,
% but may take a long time to compute.

tic

disp('Computing DOF. Please wait.')
[p_0, chisq_0,sol_0] = compute_GDF(g2,g2_error,X0s,prob,Q,N,M,L,T,w,n_perts,D,start_lm);

toc
%% Search for P_reject = 0.5 +- 0.1 with binary search

max_its=100; %terminate after max_its terations

tic

[search_lm,search_sols,search_chisq,search_P] = multiq_lm_bisection_search(g2,g2_error,Q,N,M,L,T,w,D,prob,custom_start_points,p_0,chisq_0,start_lm,max_its);

toc

%% Plot result from binary search

transform_s=1; %set to 1 to plot s in data units (velocity/diffusivity etc.)

result_best_lm = multiq_compile_result(sample,series,search_sols(:,end),q_value,t,fit_eval_t,s,T,fit_eval_T,w,q_powers,t_powers,g2,g2_error,t_norm,q_norm);

multiq_plot_result(result_best_lm,transform_s);

%% collect search results in struct

search_result = struct();
search_result.lm = search_lm;
search_result.sols = search_sols;
search_result.chisq = search_chisq;
search_result.P = search_P;
 
