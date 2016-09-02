function [M_out, F_out] = m_comb(M_in, F_in, p, q, Fec, f_fit_xy, f_fit_yy, m_fit_xx, m_fit_yy)

% Mating combinations with EDC, for a species with chromosomal sex determination.  From Cotton & Wedekind (2009)

% Calculate fecundity as rate per female, for each genotype:
M_in = sum(M_in);
F_in = Fec(:)'*F_in;

if all(M_in==0) || all(F_in==0)
    M_out = 0.*M_in;
    F_out = 0.*F_in;
else

a = M_in(1)/sum(M_in)*F_in(1)*m_fit_xx; % XXm XXf
b = M_in(1)/sum(M_in)*F_in(2)*f_fit_xy*m_fit_xx; % XXm XYf
c = M_in(1)/sum(M_in)*F_in(3)*f_fit_yy*m_fit_xx; % XXm YYf
d = M_in(2)/sum(M_in)*F_in(1); % XYm XXf
e = M_in(2)/sum(M_in)*F_in(2)*f_fit_xy; % XYm XYf 
f = M_in(2)/sum(M_in)*F_in(3)*f_fit_yy; % XYm YYf 
g = M_in(3)/sum(M_in)*F_in(1)*m_fit_yy; % YYm XXf 
h = M_in(3)/sum(M_in)*F_in(2)*f_fit_xy*m_fit_yy; % YYm XYf 
i = M_in(3)/sum(M_in)*F_in(3)*f_fit_yy*m_fit_yy; % YYm YYf

% Genotypes
XX = a + 0.5*b + 0.5*d + 0.25*e;
XY = 0.5*b + c + 0.5*d + 0.5*e + 0.5*f + g + 0.5*h;
YY = 0.25*e + 0.5*f + 0.5*h + i;

% Order of genotypes in each phenotype is XX, XY, YY
M_out = [q*XX, (1-p)*XY, (1-p)*YY];
F_out = [(1-q)*XX, p*XY, p*YY];

end