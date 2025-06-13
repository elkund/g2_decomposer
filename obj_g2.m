function [out, grad] = obj_g2(T,X,g2_data,g2_error,Q,N,M,L,lm,w,D)
% out is the objective function
% grad is the objective gradient (not implemented yet)

% inputs: T is the design matrix, X is the solution, g2s is the data, g2errs
% are the error bars, Q is the number of q bins, N is the number of delay
% points, L is the number of dynamical components, lm is the regularizer
% weight. w are the quadrature weights for regularizer integral.
% Delta_matrix is a tridiagonal matrix used to compute the
% second derivative.


% Regularizer
second_diff = D*X(1+2*Q:end);
REG = w*(second_diff.^2);

f_fits = f_gen(T,X,Q,N,M,L);
g2_fits = g2_gen(X,f_fits,Q);

% Output to be minimized.
devs = (g2_data-g2_fits)./g2_error;
out = sum(devs.^2,'all') + lm*REG;

% Compute objective gradient when requested by solver. Using this seems to
% often make results worse/convergence slower.
if nargout>1
   grad=zeros(2*Q+M*L,1);

   devs_sigma2 = devs./g2_error;
   %baseline gradientobj
   grad(1:Q) = -2*sum(devs_sigma2,2);
   %contrast gradient
   grad(1+Q:2*Q) = -2*sum(devs_sigma2.*f_fits.^2,2);
   %phi gradient (chi square)
   grad(1+2*Q:end) = -4*squeeze(sum(devs_sigma2.*X(1+Q:2*Q).*f_fits.*permute(T,[1,3,2]),[1,2]));
   %phi gradient (regularizer)
   grad(1+2*Q:end) = grad(1+2*Q:end) + 2*lm*Delta_matrix'*(w'.*second_diff);
end
end