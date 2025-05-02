function [T,w_all,s_all,D] = multiq_setup_arrays(q_value,t,s_range,q_powers,t_powers,M,varargin) 
    % input arguments: q value vectors q_value, delay value vector t,
    % integration limits s_range (shape: L by 2), per component q power
    % (shape: L), delay power (shape: L), no of points per component M,
    % optional s_scale, t_norm, q_norm (see below).

    % output arguments: design matrix T, quadrature weights w_all,
    % integration points s_all, second difference matrix D.

    % optional arguments
    p = inputParser;

    % The integration points s can be linearly or logarithmically spaced.
    % log is default.
    default_scale = 'log';

    % By default we treat t and q in units of the largest values in the
    % data set. this way, s=1 represent a decay to 1/e at t=max(t),
    % q=max(q).
    default_t_norm = max(t);
    default_q_norm = max(q_value);

    valid_scales = {'log','lin'};
    check_scale = @(x) any(validatestring(x,valid_scales));

    for arg={'q_value','t','s_range','q_powers','t_powers','M'}
        addRequired(p,arg{1},@isnumeric);
    end
    addParameter(p,'s_scale',default_scale,check_scale)
    addParameter(p,'t_norm',default_t_norm,@isnumeric)
    addParameter(p,'q_norm',default_q_norm,@isnumeric)

    p.KeepUnmatched = true;

    parse(p,q_value,t,s_range,q_powers,t_powers,M,varargin{:})

    %disp(['Scale: ',p.Results.Scale])

    % q-dependent contrast and bg
    L=length(q_powers);
    Q=length(q_value);
    N=length(t);

    w_all = zeros(1,L*M);
    s_all = zeros(1,L*M);
    T = zeros(Q,M*L,N);

    %Normalize q and t
    qc = q_value/p.Results.q_norm;
    tc = t/p.Results.t_norm;

    block_diag_cells = cell(L,1); %temporary cell array for constructing block tridiagonal matrix

    for l = 1:L
        lo = s_range(l,1);
        hi = s_range(l,2);
        switch p.Results.s_scale
            case 'lin'
                bin_edges = linspace(lo,hi,M+1);
            case 'log'
                bin_edges = logspace(log10(lo),log10(hi),M+1);
        end
        
        % We use the rectangular rule. Quadrature weights are the
        % width of each bin
        w = bin_edges(2:end) - bin_edges(1:end-1);
        %evaluation points in center of bins
        s = (bin_edges(1:end-1)+bin_edges(2:end))/2; 

        [sM,qM,tM] = meshgrid(s,qc,tc);


        T(:,1+M*(l-1):M*l,:) = w.*exp(-qM.^q_powers(l).*tM.^t_powers(l).*sM);

        w_all(1+(l-1)*M:M*l) = w;
        s_all(1+(l-1)*M:M*l) = s;

        % Construct second difference matrix.
        diff_s = diff(s)';
        
        A = 1./diff_s(1:end-1).^2;
        B = -1./diff_s(1:end-1).*(1./diff_s(1:end-1)+1./diff_s(2:end));
        C = 1./diff_s(1:end-1)./diff_s(2:end); 
    


        % Note: For lin scaled s, diff_s is constant. In that case, the
        % above could be simplified.


        block = spdiags([A B C],0:2,M-2,M);

        % zero padding to make weights and derivatives match up during
        % integration.
        block(end+1:end+2,:) = 0;
        block_diag_cells{l} = block;
    end
    
    % Block tridiagonal second difference matrx.
    D = blkdiag(block_diag_cells{:});
end