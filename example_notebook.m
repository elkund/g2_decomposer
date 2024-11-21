
year = 2024;
prop = 11019107;
sample = 'ASW_dir_CT234';
series = 32;
q_sel = 4:12;
t_sel = 0;

% Load P10 data
[qs,delays,g2s,g2errs,~,~] = p10_load(year,prop,sample,series,q_sel,t_sel);

N = length(delays); %number of observations (delay) per q bin
Q = length(qs);     %number of q bins

% Settings for plotting g2s
set(groot,'defaulttextinterpreter','latex')
xlbl = '$\tau \quad [s]$';
ylbl = '$g_2$';

%% Display data
% This is the usual g2 vs log t plot

figure();

ax1=gca;
c_map=parula(Q+1);
hold on
for q = 1:Q
    label = strjoin([compose("%5.2e",qs(q)') ' ' char(197) '^{-1}'],'');
    errorbar(delays,g2s(q,:),g2errs(q,:),'o', 'color',c_map(q,:),'DisplayName',label);
end

lgd=legend('Location','southwest');
fontsize(lgd,6,'points');
ax1.XScale="log";
xlabel(xlbl);
ylabel(ylbl);

%% MULTIQ analysis setup

% Set dynamic component power laws. For example, two components;
%   F = integral w. respect to x of
%       rho_1 e^(-x q^q_deps(1)*tau^time_deps(1))
%      +rho_2 e^(-x q^q_deps(2)*tau^time_deps(2))
% Program will solve for rho_1 and rho_2.

time_deps = [1 2];
q_deps = [2 2];

M = 100; %number of points in transformation, per component (total M*L)

L = length(q_deps);   %number of components


% Set the integration limits, i.e. s=logspace(lo,hi) or s=linspace(lo,hi)
range_los = [1E-2 1E-2];
range_his = [1E1 1E3];

% q and delay normalization. these numbers will be used to condidtion the
% problem, e^(-x(q/q_norm)^....)
q_norm = max(qs);
delay_norm = max(delays);

% Set up arrays for integral evaluation.
% Weighted kernel values T, quadrature weights w, integrand evaluation points s
[T,w,s] = setup_arrays_g2(qs,delays,range_los,range_his,q_deps,time_deps,M,'qNorm',q_norm,'tauNorm',delay_norm);

% Number of delay points in back-transformed results
fit_eval_N=100;
% Set up another T array for with more delay values. This is used to
% display fit curves that extend beyond the data points, and can be higher
% res (if fit_eval_N is large).
fit_eval_delays = logspace(log10(delays(1))-1,log10(delays(end))+1,fit_eval_N); %extended delay vector
[fit_eval_T,~,~] = setup_arrays_g2(qs,fit_eval_delays,range_los,range_his,q_deps,time_deps,M,'qNorm',q_norm,'tauNorm',delay_norm);

%% Bounds and constraints on solution

% Equality constraints.
% Fmincon will find a solution X such that Aeq*X=beq.

% Boundary conditions. Sum of lasts and firsts = 0, e.g. set slow and fast
% limits to zero.
Aeq = zeros(2,2*Q+M*L);
for l=1:L
    Aeq(2,1+2*Q+M*(l-1)) = 1;
    Aeq(2,2*Q+M*l) = 1;
end

%Constrain integral of |F| to 1 
Aeq(1,1+2*Q:end) = w;

beq = [1 0]'; %total pop. 1, boundaries = 0

% Upper and lower bounds.
% fmincon will find a solution such that lb =< X =< ub

% Bounds on component values
lb = zeros(2*Q+M*L,1); %Minimum is zero.
ub = ones(2*Q+M*L,1);
ub(2*Q+1:end) = 1./w; %every population element contriutes maximum area 1 to transformed function

% Bounds on baselines
lb(1:Q) = -0.01;
ub(1:Q) = 0.01;

% bounds on constrast
lb(Q+1:2*Q) = 0.01;
ub(Q+1:2*Q) = 0.125; %Max contrast for asw from Aerogel A and B

%% Initial guesses
% Generate a list of initial guesees
X0s = multiq_multistart_guess_gen(g2s,L,w);

%% Initialize parallell pool

parpool(64);

%% Optimization options

% Sqp alghorithm works well. Specify objective gradient not implemented yet.
% Scale problem seeems necessary even though we normalized q and tau
options = optimoptions(@fmincon,'Display','off','MaxFunctionEvaluations',3000*200,'MaxIterations',1e4,'SpecifyObjectiveGradient',false, 'Algorithm','sqp','ScaleProblem',true);

%% Set up multistart solve with fixed lagrange multiplier

lm = 1;% Manually set lagrange multiplier value

% set objective function
obj_fun  = @(x)obj_g2(T,w,x,g2s,g2errs,Q,N,s,M,L,lm);

% Create problem formulation
prob = createOptimProblem('fmincon','objective',obj_fun,'x0',X0s(:,1),'lb',lb,'ub',ub,'Aeq',Aeq,'beq',beq,'options',options);

% Create "multistart" object
ms = MultiStart('UseParallel',true);
custom_start_points = CustomStartPointSet(X0s');

%% Run multistart
% find the best solution given the start points in X0s

[Xsol1,fval,exitflag,output,all_solutions] = run(ms, prob,custom_start_points);

%% Show result

result_fixed_lm = multiq_compile_result(sample,series,Xsol1,qs,delays,fit_eval_delays,s,T,fit_eval_T,w,q_deps,time_deps,g2s,g2errs,delay_norm,q_norm);
multiq_plot_result(result_fixed_lm);

%% Start larger parallell pool

delete(gcp('nocreate'));
parpool(128);

%% Find "best" value for lagrange multiplier using bisection search

max_its = 30; % If no solution is found, search will terminate after max_its iterations
start_lm = 1E5; %starting value for lm.

% Warning: This takes a LONG time. Start as many parallell workers as
% possible.
lm_search_result = multiq_lm_bisection_search(g2s,g2errs,X0s,T,s,w,lb,ub,Aeq,beq,options,start_lm,max_its);
Xsol2 = lm_search_result.solutions(:,end);

%% Compile and plot result

result_best_lm = multiq_compile_result(sample,series,Xsol,qs,delays,fit_eval_delays,s,T,fit_eval_T,w,q_deps,time_deps,g2s,g2errs,delay_norm,q_norm);
multiq_plot_result(result_best_lm);