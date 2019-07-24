function [Egamma, Echi, Elambda, Egamma_so, WLgamma, WLchi, WLlambda, WLgamma_so] = algaasBandgaps(x, T)

if (nargin == 1)
    T = 300;
end;

c = 299.79245e6;
q = 160.21773e-21;
h = 662.60755e-36;

% the band gaps are in eV's
Egamma    = 1.519+1.155*x+0.37*x^2 - 5.41e-4*T^2/(T+204);
Echi      = 1.981+0.124*x+0.144*x^2 - 4.6e-4*T^2/(T+204);
Elambda   = 1.815+0.69*x - 6.04e-4*T^2/(T+204);
Egamma_so = 0.34-0.04*x + Egamma;

% the band edge wavelengths are in nm's
WLgamma     = c/(Egamma*q/h)*1e9;
WLchi       = c/(Echi*q/h)*1e9;
WLlambda    = c/(Elambda*q/h)*1e9;
WLgamma_so  = c/(Egamma_so*q/h)*1e9;
end