function f = fermidirac(E,Ef,T)

loadconstants;
kbeV = kb/q; % Boltzmann constant divided by eV
f = (1+exp((E-Ef)/(kbeV*T))).^(-1);

end
