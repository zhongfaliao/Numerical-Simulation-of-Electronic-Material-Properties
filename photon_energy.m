loadconstants;
wl = linspace(800,1000,201)*1e-9; %wavelength of interest in nm's
Energy = h*c./wl/(q); %corresponding photon energy in eV's

graphics = plot(wl,Energy);
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('\lambda (m)','FontSize',11);
ylabel('Energy (eV)','FontSize',11);
