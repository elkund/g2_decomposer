%% Load data and settings

%%%%%%%%
% PRE2021 313K
%%%%%%%%
load("example_data/PRE2021_313K.mat")
load("example_settings/PRE2021_313K_KWW.mat")

% For two components, set up like this:
% rate_range = [1E0,1E3;5E0,1E5]; %bounds for decay rates
% stretch_range = [0,1;1,2];      %bounds for stretch parameter
% q_powers = [2,1]; %ballistic, diffusive. We don't need to give the t
% power, it is fixed by the stretch parameter fit.

N = length(t);
Q = length(q_value);

%% Setup initial guess X0 and fit parameter bounds, lb, ub

L = length(q_powers);

% length of solution vector X = (baseline per Q,contrast per Q, strengths, decay rates, stretch)
% There are L-1 component strengths. The first is fixed to
% 1-sum(other amplitudes), as a normalizing condition.
X_size = 2*Q + L-1 + 2*L; 
X0 = zeros(X_size,1);
lb = zeros(X_size,1);
ub = ones(X_size,1);

%baselines
base_guesses = (baseline_ub+baseline_lb)/2; %initial guess for baselines
X0(1:Q) = base_guesses;
lb(1:Q) = baseline_lb; %lower bound for baselines
ub(1:Q) = baseline_ub;  %upper bound for baselines

%Contrast
contr_guesses = g2(1,1)-base_guesses-1; %contrast inital guess
X0(1+Q:2*Q) = contr_guesses;
lb(1+Q:2*Q) = contrast_lb; %contrast lower bound
ub(1+Q:2*Q) = contrast_ub; % Max contrast

%component strengths
A_slice = 1+2*Q:2*Q+L-1;
X0(A_slice) = 1/L; %iniitial guess for strengths
lb(A_slice) = 0;   %lower bound 0
ub(A_slice) = 1;   %upper bound 1

% rate, stretch parameters, X0 lb ub
G_slice = 1+2*Q+L-1:2*Q+L-1+L;
lb(G_slice) = rate_range(:,1);
ub(G_slice) = rate_range(:,2);
lb(G_slice+L) = stretch_range(:,1);
ub(G_slice+L) = stretch_range(:,2);
Gg_slice = G_slice(1):G_slice(end)+L;
X0(Gg_slice) = (lb(Gg_slice)+ub(Gg_slice))/2;

% constrain sum of strengths to =< 1.
% This is necessary and sufficient to avoid negative solutions for the first component
Aineq=zeros(size(X0))';
Aineq(1+2*Q:2*Q+L-1)=1;
bineq=1;


%% Ojective function

% obj_fun = @(X) obj_g2_KWW(X,q_powers,g2,g2_error,q_value,t,Q,L);
q_norm=max(q_value);
t_norm=max(t);
g2_fun = @(X) g2_gen(X,f_gen_KWW([1-sum(X(A_slice));X(A_slice)],X(G_slice),X(G_slice+L),q_powers,q_value/q_norm,t/t_norm),Q);
KWW_obj_fun = @(X) RSS(g2_fun(X),g2,g2_error);

%% Solve

ms_opts = optimoptions(@fmincon,'SpecifyObjectiveGradient',false,'Algorithm','sqp','MaxFunctionEvaluations',3000*200,'MaxIterations',1e4);
[sol,chisq] = fmincon(createOptimProblem('fmincon','objective',KWW_obj_fun,'x0',X0,'lb',lb,'ub',ub,'Aineq',Aineq,'bineq',bineq,'options',ms_opts));

%% Compile result and plot

KWW_result = KWW_compile_result(sample,series,sol,q_value,t,q_powers,q_norm,t_norm,g2,g2_error);

% Computing the pollard distriution involves a numerical integral. This may
% produce some warnings, but usually works.

transform_s=0;
multiq_plot_result(KWW_result,transform_s);