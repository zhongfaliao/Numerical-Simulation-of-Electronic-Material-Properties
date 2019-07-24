function delta_Eg = bandgap_shrink(x,P,N,T,E)

% reference IEEE Journal of Quantum Electronics, Vol 26, No 1, P113,
% (1990), Equation (17)
% Parameter declaraion:
% x aluminum/gallium concentration in algaas/ingaas, irrelevant for InP
% N doping
% P doping
% T absolute temperature
% E photon energy, in eV, in the equation/express, will be converted to J

% note bandgas shrinkage effects are different for n-type, p-type, and
% carrier-injection cases, they have different \kappa values
% p-type kappa = 0.11, then we should use P
% n-type kappa = 0.125, then we should use N
% carrier kappa = 0.14, use N or P, because N=P for charge neutrality
kappa = 0.125;
% however, they all share the critical carrier concentration \chi_cr

loadconstants;
props = algaas_elec_prop(x, 0, 0, T); % for AlGaAs
% props = ingaas_elec_prop(x, 0, 0, T); % for InGaAs
% props = inp_elec_prop(T); % for InP

delta_Eg = -kappa/props.eps_r_static.*((N/(props.chi_cr)-1).^(1/3)) .* heaviside(-(1-N/(props.chi_cr))); % the last number to make plot smooth

% comments: delta_Eg is negative, therefore "shrinkage"
end