function props = inp_elec_prop(T)

% Reference: IEEE Journal of Quantum Electronics, Vol 26, No 1, P113,
% (1990)

loadconstants;
kbeV = kb/q; % Boltzmann constant divided by eV

props.eps_r_static = 12.4; % Dielectric constant (static)
props.eps_r_highf = 9.61; % Dielectric constant (high frequency)

props.meff_e = 0.075; % effective electron mass
props.meff_h = 0.56; % effective (heavy) hole mass
props.meff_lh = 0.12; % effective (light) hole mass
% props.meff_cd = 0.85-0.14*x; % Density-of-states electron mass (valid when x>0.41 for AlGaAs, this parameter appears)
props.meff_vd = (props.meff_h^(3/2)+props.meff_lh^(3/2))^(2/3); % Density-of-states hole mass
props.mreduced_eh = 1/(1/props.meff_e+1/props.meff_h);
props.mreduced_elh = 1/(1/props.meff_e+1/props.meff_lh);

props.Egamma = 1.421 - 4.9e-4*T^2/(T+327);
props.Echi = 0.96 - 3.7e-4*T ;

% props.affinity = (4.9-0.83*x); % Electron affinity is very temperature insensitive
props.Eg = props.Egamma;
% props.mu_e = (40-80.7*x+49.2*x^2)*1e3; % electron Hall mobility

props.dos_c = 4.82e15*(props.meff_e*T).^1.5; % conduction band effective density of states 

% props.mu_h = 350; % hole Hall mobility
props.dos_v = 4.82e15*(props.meff_vd*T).^1.5; % valence band effective density of states, originally used mass: meff_h
props.n_i = sqrt(props.dos_c*props.dos_v)*exp(-props.Eg./(2*kbeV*T)); % Intrinsic Carrier Concentration
props.Ei = props.Eg+kbeV*T*log(props.n_i./props.dos_c);
% props.workfunc = props.Eg-props.Ei+props.affinity;
props.chi_cr = 1.6e24*(props.meff_e/(1.4*props.eps_r_static))^3; % critical electron density

% props.C = 3.38*pi/4*q^2*props.Eg/(sqrt(props.eps_r_static)*eps0*c*h^2*m0)*(2*props.mreduced_eh*m0)^(3/2)*sqrt(q);
props.C = pi/4*q^2*props.Eg/(sqrt(props.eps_r_static)*eps0*c*h^2*m0)*(2*props.mreduced_eh*m0)^(3/2)*4.40/3.05; % this is directly usable in equation 
% props.C = 7.39e5*(2*props.mreduced_eh)^(3/2)/sqrt(q);
% props.C = 6.5e5*(2*props.mreduced_eh)^(3/2)*1e8; % modified to to fit with paper.
props.C_hh = props.C*(props.mreduced_eh^(3/2)/(props.mreduced_eh^(3/2)+props.mreduced_elh^(3/2)));
props.C_lh = props.C*(props.mreduced_elh^(3/2)/(props.mreduced_eh^(3/2)+props.mreduced_elh^(3/2)));