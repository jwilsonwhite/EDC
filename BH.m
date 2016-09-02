function surv = BH(a,b,S)

% Beverton-Holt function.  Gives survival, not R

surv = a./(1 + (a/b).*S);

