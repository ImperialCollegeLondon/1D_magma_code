k=20;
smin = @(a,b) -(1/k)*log(exp(-k*a) + exp(-k*b));

smin(1500/k,800/k)