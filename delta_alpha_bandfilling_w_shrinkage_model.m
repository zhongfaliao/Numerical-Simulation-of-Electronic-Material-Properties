function delta_alpha = delta_alpha_bandfilling_w_shrinkage_model(x,Na,Nd,T,E)
% This function calculates the differential absorption coefficient in the
% reference: IEEE Journal of Quantum Electronics, Vol 26, No 1, P113, 1990
% specifically, equation (11) in the paper
% 
% Zhongfa Liao
% University of Toronto
% August 23, 2015
% 
% Declaration of parameters:
% x = Aluminum concentration in AlGaAs or Gallium concentration in GaInAs
% Na = p-doping concentration in cm-3
% Nd = n-doping concentration in cm-3
% T = temperature in Kelvin
% E = photon energy in eV

loadconstants;
props = algaas_elec_prop(x, 0, 0, T); % for AlGaAs
% props = ingaas_elec_prop(x, 0, 0, T); % for InGaAs
% props = inp_elec_prop(T); % for InP

%% Incorporate bandgap shrinkage effects here

% This determines whether to incorporate the bandgap shrinkage.

delta_Eg = bandgap_shrink(x,Na,Nd,T,E);
props.Eg = props.Eg+delta_Eg;

%% Equation (6a) (6b), calculate photon energy above/below conduction/valence band edges for electron - light/heavy holes
E_ah = (props.Eg-E)*(props.meff_e/(props.meff_e+props.meff_h))-props.Eg;
E_al = (props.Eg-E)*(props.meff_e/(props.meff_e+props.meff_lh))-props.Eg;
E_bh = (E-props.Eg)*(props.meff_h/(props.meff_e+props.meff_h));
E_bl = (E-props.Eg)*(props.meff_lh/(props.meff_e+props.meff_lh));

%% Equation (8a) (8b) y for n-type doping, z for p-type doping

y = Nd;
z = Na;

Efc = (log(y/props.dos_c)+y/props.dos_c.*(64+0.05524.*y/props.dos_c.*(64+(y/props.dos_c).^(1/2))).^(-1/4))*kb*T/q;
Efv = (-(log(z/props.dos_v)+z/props.dos_v.*(64+0.05524.*z/props.dos_v.*(64+(z/props.dos_v).^(1/2))).^(-1/4)))*kb*T/q-props.Eg;

%% Equation (7a) (7b), a function fermidirac() is defined
fc_ebh = fermidirac(E_bh,Efc,T);
fc_ebl = fermidirac(E_bl,Efc,T);
fv_eah = fermidirac(E_ah,Efv,T);
fv_eal = fermidirac(E_al,Efv,T);

%% Equation (11)
if (E<=props.Eg)
    delta_alpha = 0;
else
    delta_alpha = (props.C_hh/E*sqrt(E-props.Eg).*(fv_eah-fc_ebh-1) + props.C_lh/E*sqrt(E-props.Eg).*(fv_eal-fc_ebl-1))*sqrt(q);  
%     delta_alpha = (props.C_hh/E*sqrt(E-props.Eg).*(1) + props.C_lh/E*sqrt(E-props.Eg).*(1))*sqrt(q);
%     Energies in the above equation are in eV's, the factors are taken
%     care of by props.C_hh, lh.
%     the last factor is to calibrate J to eV conversion.
end
end