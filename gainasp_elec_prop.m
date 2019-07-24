function props = gainasp_elec_prop(x, Na, Nd, T)

% x is the concentration of Gallium in the compound semicouductor
% y is the concentration of Arsenic in the compound semicouductor
% main reference http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/index.html

loadconstants;
kbeV = kb/q; % Boltzmann constant divided by eV
props.eps_r_static = 12.5-1.4*x; % Dielectric constant (static)
props.eps_r_highf = 9.61-0.5*x; % Dielectric constant (high frequency)
props.a = 5.8687-0.4182*x;  % Lattice constant.

% props.meff_e = 0.023+0.037*x+0.003*x^2; % effective electron mass
props.meff_h = 0.6+0.19*x; % effective (heavy) hole mass
props.meff_lh = 0.09+0.05*x; % effective (light) hole mass
props.meff_cd = 0.63+0.13*x; % Density-of-states electron mass (valid when x>0.74 for GaxIn(1-x)AsyP(1-y), this parameter appears)
% props.meff_vd = (props.meff_h^(3/2)+props.meff_lh^(3/2))^(2/3); % Density-of-states effective hole mass
% props.mreduced_eh = 1/(1/props.meff_e+1/props.meff_h);
% props.mreduced_elh = 1/(1/props.meff_e+1/props.meff_lh);

% props.Egamma = 0.4105+0.6337*x+0.475*x^2 -(6*x^2- 8.6*x +5.2)*1e-4*T^2/(337*x^2-455*x+196); % on GaAs substrate, Temperature
% props.Egamma = 0.42 + 0.625*x -(5.8/(T+300)-4.19/(T+271))*1e-4*T^2*x-4.19e-4*T^2/(T+271) +0.475*x^2;
% props.Echi = 1.37-0.63*x+1.16*x^2; % 300k, temperature dependence info not available

% props.affinity = 4.38-0.58*x; % Electron affinity is very temperature insensitive
% props.Eg = props.Egamma;
% props.mu_e = (40-80.7*x+49.2*x^2)*1e3; % electron Hall mobility

% props.dos_c = 4.82e15*(props.meff_e*T).^1.5; % conduction band effective density of states 

% props.mu_h = 350; % hole Hall mobility
% props.dos_v = 4.82e15*(props.meff_vd*T).^1.5; % valence band effective density of states, originally used mass: meff_h
% props.n_i = sqrt(props.dos_c*props.dos_v)*exp(-props.Eg./(2*kbeV*T)); % Intrinsic Carrier Concentration
% props.Ei = props.Eg+kbeV*T*log(props.n_i./props.dos_c);
% props.workfunc = props.Eg-props.Ei+props.affinity;
% props.chi_cr = 1.6e24*(props.meff_e/(1.4*props.eps_r_static))^3; % critical electron density

% props.C = 3.38*pi/4*q^2*props.Eg/(sqrt(props.eps_r_static)*eps0*c*h^2*m0)*(2*props.mreduced_eh*m0)^(3/2)*sqrt(q);
props.C = pi/4*q^2*props.Eg/(sqrt(props.eps_r_static)*eps0*c*h^2*m0)*(2*props.mreduced_eh*m0)^(3/2);
% props.C = 7.39e5*(2*props.mreduced_eh)^(3/2)/sqrt(q);
% props.C = 6.5e5*(2*props.mreduced_eh)^(3/2)*1e8; % modified to to fit with paper.
props.C_hh = props.C*(props.mreduced_eh^(3/2)/(props.mreduced_eh^(3/2)+props.mreduced_elh^(3/2)));
props.C_lh = props.C*(props.mreduced_elh^(3/2)/(props.mreduced_eh^(3/2)+props.mreduced_elh^(3/2)));