function [A,w_all,s_all] = setup_arrays_g2(qs,tau,range_los,range_his,q_deps,time_deps,M,varargin) 

    % optional arguments
    p = inputParser;

    defaultScale = 'log';
    defaulttauNorm = max(tau);
    defaultqNorm = max(qs);

    validScales = {'log','lin'};
    checkScale = @(x) any(validatestring(x,validScales));

    for arg={'qs','tau','range_los','range_his','q_deps','time_deps','M'}
        addRequired(p,arg{1},@isnumeric);
    end
    addParameter(p,'Scale',defaultScale,checkScale)
    addParameter(p,'tauNorm',defaulttauNorm,@isnumeric)
    addParameter(p,'qNorm',defaultqNorm,@isnumeric)

    p.KeepUnmatched = true;

    parse(p,qs,tau,range_los,range_his,q_deps,time_deps,M,varargin{:})

    %disp(['Scale: ',p.Results.Scale])

    % q-dependent contrast and bg
    L=length(q_deps);
    Q=length(qs);
    N=length(tau);

    w_all = zeros(1,L*M);
    s_all = zeros(L,M);
    A = zeros(Q,M*L,N);

    %Normalize q and tau
    qsc = qs/p.Results.qNorm;
    tauc = tau/p.Results.tauNorm;

    for l = 1:L
        switch p.Results.Scale
            case 'lin'
                lo = range_los(l);
                hi = range_his(l);
                bin_edges = linspace(lo,hi,M+1);
            case 'log'
                lo = log10(range_los(l));
                hi = log10(range_his(l));
                bin_edges = logspace(lo,hi,M+1);
        end
        
        % We use the rectangular rule. Quadrature weights are the
        % width of each bin
        w = bin_edges(2:end) - bin_edges(1:end-1);
        %evaluation points in center of bins
        s = (bin_edges(1:end-1)+bin_edges(2:end))/2; 

        [sM,qM,tM] = meshgrid(s,qsc,tauc);

        A(:,1+M*(l-1):M*l,:) = w.*exp(-qM.^q_deps(l).*tM.^time_deps(l).*sM);

        w_all(1+(l-1)*M:M*l) = w;
        s_all(l,:) = s;
    end
end