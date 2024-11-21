function [out, grad] = obj_g2(T,w,X,g2s,g2errs,Q,N,s,M,L,lm)
% out is the objective functoin to be minimized.
% grad is the objective gradient (not implemented yet)

%regularizer
REG = 0;
if lm ~=0
    %conditioning = median(s,'all');
    for l=1:L
        %some of this stuff doesn't need to be computed every time.
         x=s(l,:);
         w_sub = w(1+M*(l-1):M*l);
         y = X(1+2*Q+M*(l-1):2*Q+M*l);

         dydx = diff(y)./diff(x)';
         x1 = (x(1:end-1)+x(2:end))/2;
         w2 = diff(x1);
         d2ydx2 = diff(dydx)./w2';

         %Computing the 2nd with central difference leads
         %to possible jagged solutions; e.g. y=[1 0 1 0 1 0] gives small
         % curvature at every point. Changed on 10.10.2024
         %r = DGradient(DGradient(y,x,[],'1stOrder'),x,[],'1stOrder');%diff(diff(g(1+2*K+M*(l-1):2*K+M*l));

         norm = w_sub*y;
         if norm ~= 0
             reg = w2*(d2ydx2.^2);
             REG = REG + reg/norm;
         end
    end
else
    REG=0;
end

g2fits=g2_gen(X,f_gen(T,X,Q,N,M,L),Q);

% Output to be minimized.
out = MSD(g2fits,g2s,g2errs)+lm*REG;

% TODO: compute objective gradient
%i nargin>1 ??
%   grad=?
end