function delta_n = delta_n_free_carrier_absorption_model(x,P,N,T,E)

% reference IEEE Journal of Quantum Electronics, Vol 26, No 1, P113,
% (1990), Equation (21)
% Parameter declaraion:
% x aluminum/gallium concentration in algaas/ingaas, irrelevant for InP
% N doping
% P doping
% T absolute temperature
% E photon energy, in eV, in the equation/express, will be converted to J

loadconstants;
props = algaas_elec_prop(x, 0, 0, T); % for AlGaAs
% props = ingaas_elec_prop(x, 0, 0, T); % for InGaAs
% props = inp_elec_prop(T); % for InP

delta_n = -6.9e-22/sqrt(props.eps_r_static)/E^2*(N/props.meff_e+P*(props.meff_h^0.5+props.meff_lh^0.5)/(props.meff_h^1.5+props.meff_lh^1.5)); %/(q^2*m0);

% comment: the value of delta_n is always negative.

end