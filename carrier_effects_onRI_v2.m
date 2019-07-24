% Carrier Induced Change in Refractive Index of InP, GaAs and InGaAsP
% Zhongfa Liao
% University of Toronto
% Aug 18th, 2015
% This code calculates III-V semiconduct refractive index change due to
% carrier injections and doping, this version (v2) compared to the last
% version, has the ability to calculation doping, and carrier injection
% separately.
% Reference: IEEE Journal of Quantum Electronics, Vol 26, No 1, P113,
% (1990)

%% Initialize
clear all;clc;close all;

%% What material to simulate?
% In this step, we need to pick a material to simulate.

x = 0.0;
% Aluminum/Gallium content in AlGaAs/InGaAs Note, we also need to change
% corresponding functions to get material properties in Step (1), Step (2)
% and Step (4).
% x is irrelevant when we simulate InP.

T = 293;
% temperature in Kelvin, therefore T = 293K corresponds to T = 20C.
E = 1.0:0.001:2.0;
% photon enery in eV
% Note, the photon energy range needs to be larger than the range in which
% we want to get the change of refractive index or absorption coefficient.
% This is due to the fact that we are to use Kramers-Kronig relation.

% y and z could be n-doping, p-doping, or n- and p- doping (for carrier
% injection/depletion case.)
z = cat(2, linspace(1E17,5E17,2),linspace(1E18,5E18,2));
y = 5e15; % This corresponds to un-intentional doping level.
% y = z; % For carrier injection/depletion, it should be N = P.

%% Step 1(0) Bandfilling
% This step calculates bandfilling effects with considering band gap
% shrinkage effects, which is reflected in the delta_alpha_bandfilling
% function. In delta_alpha_bandfilling, band gap shrinkage function is
% disabled.

delta_alpha_bandfilling = NaN(length(z),length(E));
% delta_alpha data for different carrier concentrations are stored in
% different rows; delta_alpha data for different photon energies
% (wavelengths) are stored in different columns.

for m = 1:length(z)
    for n = 1:length(E)
        delta_alpha_bandfilling(m,n) = delta_alpha_bandfilling_model(x,z(m),y,T,E(n));
    end
end

figure('Name','Band filling effect (w/o band gap shrinkage effect)','NumberTitle','on');
graphics = semilogy(E, -delta_alpha_bandfilling(1,:),E, -delta_alpha_bandfilling(2,:), E, -delta_alpha_bandfilling(3,:), E, -delta_alpha_bandfilling(4,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('-\Delta\alpha (cm^-^1)','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
ylim([10^1 10^5]);
legend('N = 1\times10^1^7 cm^-^3', 'N = 5\times10^1^7 cm^-^3', 'N = 1\times10^1^8 cm^-^3','N = 5\times10^1^8 cm^-^3','Location','northwest');

%% Step 1(1) Base absorption coefficient without doping/carrier injection
% Base absorption coefficient is calculated when the carrier concentration
% is at the un-intentional doping level, which is 5E15 cm^(-3).

alpha0_bandfilling = zeros(1,length(E));

for n=1:length(E)
    alpha0_bandfilling(n) = delta_alpha_bandfilling_model(x,5E15,5E15,T,E(n));
end

%% Step 1(2) Kramers-Kronig relation
% Base absorption coefficient calculated in the above step (1) is for
% simplifying the normalization process.

delta_n_bandfilling = NaN(length(z),length(E));

for k = 1:length(z)
    delta_n_bandfilling(k,:) = kkintegral(E,delta_alpha_bandfilling(k,:)-alpha0_bandfilling,E);
end

figure('Name','Band filling effect (w/o band gap shrinkage effect)','NumberTitle','on');
graphics = plot(E,delta_n_bandfilling(1,:),E,delta_n_bandfilling(2,:),E,delta_n_bandfilling(3,:),E,delta_n_bandfilling(4,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('\Deltan','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([-0.005 0.005]);
legend('N = 1\times10^1^7 cm^-^3', 'N = 5\times10^1^7 cm^-^3', 'N = 1\times10^1^8 cm^-^3','N = 5\times10^1^8 cm^-^3','Location','northwest');

%% Step 2(0) Bandgap shrinkage

delta_alpha_bandgap_shrinkage = NaN(length(z),length(E));

for m = 1:length(z)
    for n = 1:length(E)
        delta_alpha_bandgap_shrinkage(m,n) = delta_alpha_bandgap_shrinkage_model(x,z(m),y,T,E(n));
    end
end

figure('Name','Band gap shrinkage effect','NumberTitle','on');
graphics = plot(E, delta_alpha_bandgap_shrinkage(1,:), E, delta_alpha_bandgap_shrinkage(2,:), E, delta_alpha_bandgap_shrinkage(3,:), E, delta_alpha_bandgap_shrinkage(4,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('\Delta\alpha (cm^-^1)','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([10^1 10^5]);
legend('N = 1\times10^1^7 cm^-^3', 'N = 5\times10^1^7 cm^-^3', 'N = 1\times10^1^8 cm^-^3','N = 5\times10^1^8 cm^-^3','Location','northwest');

%% Step 2(1) Initial absorption without doping/injection

alpha0_bandgap_shrinkage = zeros(1,length(E));
for n=1:length(E)
    alpha0_bandgap_shrinkage(n) = delta_alpha_bandgap_shrinkage_model(x,5E15,5E15,T,E(n));
end

%% Step 2(2) Kramers-Kronig relation

delta_n_bandgap_shrinkage = NaN(length(z),length(E));

for k = 1:length(z)
    delta_n_bandgap_shrinkage(k,:) = kkintegral(E,delta_alpha_bandgap_shrinkage(k,:)-alpha0_bandgap_shrinkage,E);
end

figure('Name','Band gap shrinkage effect','NumberTitle','on');
graphics = plot(E,delta_n_bandgap_shrinkage(1,:),E,delta_n_bandgap_shrinkage(2,:),E,delta_n_bandgap_shrinkage(3,:),E,delta_n_bandgap_shrinkage(4,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('\Deltan','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([-0.101 0.101]);
legend('N = 1\times10^1^7 cm^-^3', 'N = 5\times10^1^7 cm^-^3', 'N = 1\times10^1^8 cm^-^3','N = 5\times10^1^8 cm^-^3','Location','northwest');

%% Step 3 Free carrier absorption

% This step directly calculates the refractive index change.
% Kramers-Kronig theory was not used.

delta_n_free_carrier_absorption = NaN(length(z),length(E));

for i=1:length(z)
    for j=1:length(E)
        delta_n_free_carrier_absorption(i,j) = delta_n_free_carrier_absorption_model(x,z(i),y,T,E(j));
    end
end

figure('Name','Free carrier absorption effect','NumberTitle','on');
graphics = plot(E,delta_n_free_carrier_absorption(1,:),E,delta_n_free_carrier_absorption(2,:),E,delta_n_free_carrier_absorption(3,:),E,delta_n_free_carrier_absorption(4,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('-\Deltan','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([10^1 10^5]);
legend('N = 1\times10^1^7 cm^-^3', 'N = 5\times10^1^7 cm^-^3', 'N = 1\times10^1^8 cm^-^3','N = 5\times10^1^8 cm^-^3','Location','northwest');

%% Step 4(0) Re-calculate bandfilling effect when band gap shrinkage is considered.

delta_alpha_bandfilling_with_shrinkage = NaN(length(z),length(E));

for m = 1:length(z)
    for n = 1:length(E)
        delta_alpha_bandfilling_with_shrinkage(m,n) = delta_alpha_bandfilling_w_shrinkage_model(x,z(m),y,T,E(n));
    end
end

figure('Name','Band filling effect (with band gap shrinkage effect)','NumberTitle','on');
graphics = semilogy(E, -delta_alpha_bandfilling_with_shrinkage(1,:),E, -delta_alpha_bandfilling_with_shrinkage(2,:), E, -delta_alpha_bandfilling_with_shrinkage(3,:), E, -delta_alpha_bandfilling_with_shrinkage(4,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('-\Delta\alpha (cm^-^1)','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([10^1 10^5]);
legend('N = 1\times10^1^7 cm^-^3', 'N = 5\times10^1^7 cm^-^3', 'N = 1\times10^1^8 cm^-^3','N = 5\times10^1^8 cm^-^3','Location','northwest');

%% Step 4(1) Initial absorption without doping/injection

alpha0_bandfilling_with_shrinkage = zeros(1,length(E));

for n=1:length(E)
    alpha0_bandfilling_with_shrinkage(n) = delta_alpha_bandfilling_w_shrinkage_model(x,5E15,5E15,T,E(n));
end

%% Step 4(2) Kramers-Kronig relation

delta_n_bandfilling_with_shrinkage = NaN(length(z),length(E));

for k = 1:length(z)
    delta_n_bandfilling_with_shrinkage(k,:) = kkintegral(E,delta_alpha_bandfilling_with_shrinkage(k,:)-alpha0_bandfilling_with_shrinkage,E);
end

figure('Name','Band filling effect (with band gap shrinkage effect)','NumberTitle','on');
graphics = plot(E,delta_n_bandfilling_with_shrinkage(1,:),E,delta_n_bandfilling_with_shrinkage(2,:),E,delta_n_bandfilling_with_shrinkage(3,:),E,delta_n_bandfilling_with_shrinkage(4,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('\Deltan','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([-0.005 0.005]);
legend('N = 1\times10^1^7 cm^-^3', 'N = 5\times10^1^7 cm^-^3', 'N = 1\times10^1^8 cm^-^3','N = 5\times10^1^8 cm^-^3','Location','northwest');

%% Collective effect on the change of RI of Step (4), Step (2) and Step (3).

% Collective effect of band filling (with band gap shrinkage), band gap
% shrinkage, and free carrier absorption.
delta_n_collective = delta_n_bandfilling_with_shrinkage+ delta_n_bandgap_shrinkage  + delta_n_free_carrier_absorption;

% delta_n_collective = delta_n_bandfilling + delta_n_bandgap_shrinkage  + delta_n_free_carrier_absorption;
% Note, when calculating carrier injection/depletion induced change in the
% refractive index, band gap shrinkage effect should not be considered in
% the calculation of band filling effect. Otherwise there would be a
% dicrepancy/discruption in the collective change of delta n.

figure('Name','Collective effect on the change of RI','NumberTitle','on');
graphics = plot(E,delta_n_collective(1,:),E,delta_n_collective(2,:),E,delta_n_collective(3,:),E,delta_n_collective(4,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('\Deltan','FontSize',15);
xlim([1.2 1.501]);
% xlim([0.8 2.0]);
% xlim([min(E) max(E)]);
% ylim([-0.1 0.1]);
legend('N = 1\times10^1^7 cm^-^3', 'N = 5\times10^1^7 cm^-^3', 'N = 1\times10^1^8 cm^-^3','N = 5\times10^1^8 cm^-^3','Location','northwest');

%% Plot different effects (band filling, band gap shrinkage and free carrier absorption) in one figure.

figure('Name','Plot different effects in one figure','NumberTitle','on');
graphics = plot(E,delta_n_bandfilling(3,:),E,delta_n_bandgap_shrinkage(3,:),E,delta_n_free_carrier_absorption(3,:),E,delta_n_collective(3,:));
% Note, (3,:) in the above plot means we have chosen the third carrier
% concentration, which is 1E18 in our case. We can change to other
% concentrations.
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('\Deltan','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([-0.0601 0.0601]);
legend('Band filling', 'Band gap shrinkage', 'Free carrier absorption','Collective','Location','northwest');

%% Plot different effects on the change of absorption coefficient in one figure.

collective_delta_alpha = delta_alpha_bandfilling_with_shrinkage+delta_alpha_bandgap_shrinkage;
% Collective effect on the change of absorption of Step (4) and Step (2).
% Because Step (3) does not allude to aborption, it is irrelevant.

% Note, in this part, we always band filling effect that is calculated with
% band gap shrinkage effect considered.

figure('Name','Plot different effects on the change of absorption coefficient in one figure','NumberTitle','on');

graphics = plot(E,delta_alpha_bandfilling_with_shrinkage(3,:),E,delta_alpha_bandgap_shrinkage(3,:),E,collective_delta_alpha(3,:));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('\Delta\alpha (cm^-^1)','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([-4E4 3E4]);
legend('Band filling', 'Band gap shrinkage','Collective','Location','southwest');

%% Plot different effects on the change of imaginary part of RI in one figure.

collective_delta_alpha = delta_alpha_bandfilling_with_shrinkage+delta_alpha_bandgap_shrinkage;
% Collective effect on the change of absorption of Step (4) and Step (2).
% Because Step (3) does not allude to aborption, it is irrelevant.

loadconstants;
alpha_to_img_RI = 100*h*c./(E*q^2)/4/pi;
% The above equation translates changes in absorption coefficient to
% imaginary part of RI.

figure('Name','Plot different effects on the change of imaginary part of RI','NumberTitle','on');

graphics = plot(E,abs(delta_alpha_bandfilling_with_shrinkage(2,:)./alpha_to_img_RI),E,abs(delta_alpha_bandgap_shrinkage(2,:)./alpha_to_img_RI),E,abs(collective_delta_alpha(2,:)./alpha_to_img_RI));
set(gca,'FontSize',13);
set(gca,'LineWidth',2.5);
set(graphics,'LineWidth',2.5);
xlabel('Photon energy (eV)','FontSize',15);
ylabel('\Deltan_I_M_G','FontSize',15);
xlim([1.2 1.501]);
% xlim([min(E) max(E)]);
% ylim([-0.15 0.15]);
legend('Band filling', 'Band gap shrinkage','Collective','Location','southwest');

%% End of code