function props = algaas_elec_prop(x, Na, Nd, T)

% Main reference: http://www.ioffe.rssi.ru/SVA/NSM/Semicond/AlGaAs/index.html
% Reference: IEEE Journal of Quantum Electronics, Vol 26, No 1, P113,
% (1990)
% some of the parameters are added accordingly


loadconstants;
kbeV = kb/q; % Boltzmann constant divided by eV
props.eps_r_static = 12.9-2.84*x; % Dielectric constant (static)
props.eps_r_highf = 10.89-2.73*x; % Dielectric constant (high frequency)

props.meff_e = 0.063+0.083*x; % effective electron mass
props.meff_h = 0.51+0.25*x; % effective (heavy) hole mass
props.meff_lh = 0.082+0.068*x; % effective (light) hole mass
props.meff_cd = 0.85-0.14*x; % Density-of-states electron mass (valid when x>0.41)
props.meff_vd = (props.meff_h^(3/2)+props.meff_lh^(3/2))^(2/3); % Density-of-states effective hole mass
props.mreduced_eh = 1/(1/props.meff_e+1/props.meff_h);
props.mreduced_elh = 1/(1/props.meff_e+1/props.meff_lh);

props.Egamma = 1.424+1.247*x-5.405e-4*(T^2/(T+204)-300^2/(204+300));
props.Echi = 1.9+0.125*x+0.143*x^2-5.405e-4*(T^2/(T+204)-300^2/(204+300));

if (x < 0.45)
    props.affinity = 4.07-0.7482*x+2.702e-4*T.^2/(T+204);
    props.Eg = props.Egamma; % Egamma dominates
    props.mu_e = 8000-22000*x+10000*x^2; % electron Hall mobility
else
    props.affinity = 3.594+0.3738*x-0.143*x.^2+2.702e-4*T^.2/(T-204);
    props.Eg = props.Echi; % Echi dominates
    props.mu_e = 255+1160*x-720*x^2;
end;

if (x < 0.41)
    props.dos_c = 4.82e15*(props.meff_e*T).^1.5; % conduction band effective density of states 
else
    props.dos_c = 4.82e15*(props.meff_cd*T).^1.5;
end;

props.mu_h = 370-970*x+740*x^2; % hole Hall mobility
props.dos_v = 4.82e15*(props.meff_vd*T).^1.5; % valence band effective density of states, originally used mass: meff_h
props.n_i = sqrt(props.dos_c*props.dos_v)*exp(-props.Eg./(2*kbeV*T)); % Intrinsic Carrier Concentration
props.Ei = props.Eg+kbeV*T*log(props.n_i./props.dos_c);
props.workfunc = props.Eg-props.Ei+props.affinity;
props.chi_cr = 1.6e24*(props.meff_e/(1.4*props.eps_r_static))^3; % critical electron density

% props.C = 3.38*pi/4*q^2*props.Eg/(sqrt(props.eps_r_static)*eps0*c*h^2*m0)*(2*props.mreduced_eh*m0)^(3/2)*sqrt(q);
props.C = pi/4*q^2*props.Eg/(sqrt(props.eps_r_static)*eps0*c*h^2*m0)*(2*props.mreduced_eh*m0)^(3/2);
% props.C = 7.39e5*(2*props.mreduced_eh)^(3/2)/sqrt(q);
% props.C = 6.5e5*(2*props.mreduced_eh)^(3/2)*1e8; % modified to to fit with paper.
props.C_hh = props.C*(props.mreduced_eh^(3/2)/(props.mreduced_eh^(3/2)+props.mreduced_elh^(3/2)));
props.C_lh = props.C*(props.mreduced_elh^(3/2)/(props.mreduced_eh^(3/2)+props.mreduced_elh^(3/2)));