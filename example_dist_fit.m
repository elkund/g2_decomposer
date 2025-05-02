%% Load data and settings

%%%%%%%%%%%
% PRE2021 313 K
%%%%%%%%%%%
load("example_data/PRE2021_313K.mat")
load("example_settings/PRE2021_313K_dist_fit.mat")

% %For two components, set up like this (for example):
% median_range = [1E0,1E3;5E0,1E5]; %bounds for distribution median, exp(mu)
% sigma_range = [1E-2,5E-1;1E-2,5E-1];  %bounds for sigma parameter
% q_powers = [2,2];
% t_powers = [1,2];


N = length(t);
Q = length(q_value);

%% Set up dynamical distributions

% We use Matlab distribution objects to handle the density function.

L = length(t_powers);

dist_strs = cell(1,L);
dist_strs(:) = {'LogNormal'}; %Choose what pdfs to use. Here we set all to log normal

dists = cell(1,L);
dist_Nparams = zeros(1,L);

for l=1:L
    dists{l} = makedist(dist_strs{l}); %Hold the distribution objects
    dist_Nparams(l) = dists{l}.NumParameters; %the number of distribution parameters per component
end

%% Setup initial guess X0 and fit parameter bounds, lb, ub

% length of solution vector X = (baseline per Q,contrast per Q, strengths, dist params)
% There are L-1 component amplitudes/strengths. The first is fixed to
% 1-sum(other amplitudes), as a normalizing condition.
X_size = 2*Q + L-1 + sum(dist_Nparams); 
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
X0(1+2*Q:2*Q+L-1) = 1/L; %iniitial guess for strengths
lb(1+2*Q:2*Q+L-1) = 0;   %lower bound 0
ub(1+2*Q:2*Q+L-1) = 1;   %upper bound 1

%Dist parameters, X0 lb ub
for l = 1:L
    params_start_ind = 1+2*Q+L-1+sum(dist_Nparams(1:l-1)); %index of first dist. parameter for component l
    
    min_mu = log(median_range(l,1));
    max_mu = log(median_range(l,2));
    min_sigma = sigma_range(l,1);
    max_sigma = sigma_range(l,2);

    X0(params_start_ind) = (max_mu+min_mu)/2; % Init guess in middle of range
    X0(params_start_ind+1) = (max_sigma+min_sigma)/2; %initial guess for sigma (shape)
    
    lb(params_start_ind) = min_mu;      %-1*log(10); %lower bound on mu
    lb(params_start_ind+1) = min_sigma; %lower bound on sigma
    ub(params_start_ind) = max_mu;      %upper bound on mu
    ub(params_start_ind+1) = max_sigma; %upper bound on sigma
end

% constrain sum of strengths to =< 1.
% This is necessary and sufficient to avoid negative solutions for the first component
Aineq=zeros(size(X0))';
Aineq(1+2*Q:2*Q+L-1)=1;
bineq=1;

dist_Nx = 30; %number of points in integration (0 to 99th percentile)

% Conditions the problem with normalized q and t
q_norm = max(q_value);
t_norm = max(t);
%% Ojective function

obj_fun = @(X) obj_g2_dist_fit(X,q_powers,t_powers,dist_strs,dist_Nx,dist_Nparams,g2,g2_error,q_value/q_norm,t/t_norm,Q,N,L);

%% Solve

ms_opts = optimoptions(@fmincon,'Algorithm','sqp','MaxFunctionEvaluations',3000*200,'MaxIterations',1e4);
[dist_sol,chisq] = fmincon(createOptimProblem('fmincon','objective',obj_fun,'x0',X0,'lb',lb,'ub',ub,'Aineq',Aineq,'bineq',bineq,'options',ms_opts));

%% Compile result as a struct and plot

dist_result = dist_fit_compile_result(sample,series,dist_sol,q_value,t,q_powers,t_powers,q_norm,dist_strs,dist_Nx,dist_Nparams,lb,ub,t_norm,g2,g2_error);

transform_s=1;
multiq_plot_result(dist_result,transform_s);
