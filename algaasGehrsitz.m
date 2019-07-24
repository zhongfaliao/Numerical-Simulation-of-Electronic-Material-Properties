% Gehrsitz model AlGaAs Refractive Index Calculator
%
% Sean Wagner
% October 4, 2005
% 
% $Id: algaasGehrsitz.m,v 1.1.1.1 2010-01-14 14:34:09 wagner Exp $
%
% Based on mathematical formulae (with corrections) given in:
% S. Gehrsitz et al., "The refractive index of AlGaAs below the band gap:
% Accurate determination and empirical modeling," Journal of Applied Physics,
% vol 87, pp. 7825-7837, 2000.
%
% Function: algaasGehrsitz
%   Calculates the refractive index for a given wavelength and composition
%   fraction.
%
%   Parameters:
%       x       - fractional composition of Al (composition of Ga is 1-x)
%       lambda   - free-space wavelength (micrometers)
%				T				- temperature (Kelvin)
%
%   Returns:
%       n       - Refractive index
%

function n = algaasGehrsitz(x, lambda, T)

% Constants
h 	= 4.135669240e-15;			% Plank's constant (eV*s)
c 	= 2.997924580e8;				% Speed of light (m)
kB 	= 8.617385692e-5;				% Boltzmann constant (meV/K)

% Conversions
eV2invum = 1/h/c/1e6;

% Precalculations
E = 1./lambda;

%
% R(x,E)
%
% coefficients
C2GaAs	= 1.55e-3;
E2GaAs	= sqrt(0.724*1e-3);
C2AlAs	= 2.61e-3;
E2AlAs	= sqrt(1.331*1e-3);
% calculations
R = (1-x).*C2GaAs./(E2GaAs^2 - E.^2) + x.*C2AlAs./(E2AlAs^2 - E.^2);

%
% A(x,lambda)
%
% coefficients
c11 = -16.159;
c12 = 43.511;
c13 = -71.317;
c14 = 57.535;
c15 = -17.451;
a0	= 5.9613;			% FIT 2 param for GaAs
a1	= 7.178e-4;		% FIT 2 param for GaAs
a2	= -0.953e-6;	% FIT 2 param for GaAs
% calculations
A0 	= a0 + a1*T + a2*T.^2;			% A(0, T), i.e. GaAs
c10	= A0;
A		= c10 + c11*x + c12*x.^2 + c13*x.^3 + c14*x.^4 + c15*x.^5;

%
% C0
%
% coefficients
c40	= 50.535;
c41 = -150.7;
c42 = -62.209;
c43 = 797.16;
c44 = -1125;
c45 = 503.79;
% calculations
invC0 = c40 + c41*x + c42*x.^2 + c43*x.^3 + c44*x.^4 + c45*x.^5;
C0 = 1./invC0;

%
% E0
%
% coefficients
Egamma0 = 1.5192;						% (eV)
EDeb		=	15.9e-3;					% (eV)
ET0			= 33.6e-3;					% (eV)
S				= 1.8;
ST0			= 1.1;
c51			= 1.1308;
c52			= 0.1436;
% calculations
EgammaGaAs = Egamma0+S*EDeb*(1-coth(EDeb/2/kB./T))+ST0*ET0*(1-coth(ET0/2/kB./T));
c50 = EgammaGaAs*eV2invum;
E0 = c50 + c51*x + c52*x.^2;
%Egamma = EgammaGaAs + 1.36*x + 0.22*x^2;
%E0 = Egamma;

%
% C1
%
% coefficients
c20	= 21.5647;
c21	= 113.74;
c22 = -122.5;
c23 = 108.401;
c24 = -47.318;
% calculations
C1 = c20 + c21*x + c22*x.^2 + c23*x.^3 + c24*x.^4;

%
% E1
%
% coefficients
e0	= 4.7171;				% FIT 2 for GaAs
e1	= -3.237e-4;		% FIT 2 for GaAs
e2	= -1.358e-6;		% FIT 2 for GaAs
c31 = 11.006;
c32 = -3.08;
% calculations
E10sqrd = e0 + e1*T + e2*T.^2;
c30 = E10sqrd;
E1sqrd = c30 + c31*x + c32*x.^2;

%
% n
%
nsqrd = A + C0./(E0.^2-E.^2) + C1./(E1sqrd-E.^2) + R;
n = sqrt(nsqrd);
