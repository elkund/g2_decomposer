function out = pollard(s,g)

fun = @(x) 1/pi * sin(x.^g*sin(pi*g)).*exp(-x.^g*cos(g*pi)-s*x);
out = integral(fun,0,inf,'ArrayValued',true);

end