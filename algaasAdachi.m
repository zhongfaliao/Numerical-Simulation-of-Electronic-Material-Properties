function [n k] = algaasAdachi(x, lambda, T);

loadconstants;

Egamma    = 1.519+1.155*x+0.37*x^2 - 5.41e-4*T^2/(T+204);
Echi      = 1.981+0.124*x+0.144*x^2 - 4.6e-4*T^2/(T+204);
Elambda   = 1.815+0.69*x - 6.04e-4*T^2/(T+204);
Egamma_so = 0.34-0.04*x;

E1vb_depth = 1;
E2vb_depth = 2.80;

Apoly = [59.26 -24.59 15.76 3.239];
B1poly = [4.258 -5.814 0.3816 6.392]

params.E0 = Egamma;
params.delta0 = Egamma_so;
params.E1 = Elambda+E1vb_depth;
params.delta1 = 3.13-params.E1;
params.E2 = Echi+E2vb_depth;
if (x<0.4)
  params.EgID = Elambda;
else
  params.EgID = Echi;
end
params.A = polyval(Apoly, x);
params.B1 = polyval(B1poly, x);
params.Gamma = 0.11;
params.C = 2.39;
params.gamma = 0.146;
params.D = 24.2;
params.omegaq = 0;
alpha0 = 5.65; 

omega = 2*pi*c/(lambda*1e-6);

params.B2 = 44*(params.E1+2*params.delta1/3)/(alpha0*(params.E1+params.delta1)^2);

[n k] = dispAdachi(omega, params);
