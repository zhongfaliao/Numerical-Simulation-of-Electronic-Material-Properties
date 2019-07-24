% Zhongfa Liao
% University of Toronto
% Aug 25, 2015
% The carrier effects on the change of refractive index for n?type GaAs at
% wavelength = 1.06,1.3, and 1.55 um
% Journal of Applied Physics 67, 1497 (1990);

clear all;
clc;
close all;

loadconstants;
x = 0.0:0.01:1.0; % Aluminum concentration in AlGaAs

%% Plot energy band gaps for different x
Egamma = NaN(1, length(x));
Echi = NaN(1, length(x));
Elambda = NaN(1, length(x));
Egamma_so = NaN(1, length(x));

WLgamma = NaN(1, length(x));
WLchi = NaN(1, length(x));
WLlambda = NaN(1, length(x));
WLgamma_so = NaN(1, length(x));

for i=1:length(x)
    [Egamma(i), Echi(i), Elambda(i), Egamma_so(i), WLgamma(i), WLchi(i), WLlambda(i), WLgamma_so(i)] = algaasBandgaps(x(i));
end

figure;
graphics = plot(x, Egamma, x, Echi, x, Elambda, x, Egamma_so);
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Aluminum content ','FontSize',11);
ylabel('Band gap energy (eV)','FontSize',11);
legend('E\Gamma','E\chi','E\Lambda','E\Gammaso','Location','northwest');

%% Band filling effects on the shift of band gap
% Due to band filling, effective Fermi energy level shifts to higher value,
% or put it in another way, band gap shifts to higher value
% Therefore, Eg = Eg0 + EF

y=cat(2, linspace(1e16,9e16,9),linspace(1e17,9e17,9),linspace(1e18,9e18,9),linspace(1e19,1e20,10));

props0 = algaas_elec_prop(0.0, 0, 0, 300);
props1 = algaas_elec_prop(0.1, 0, 0, 300);
props2 = algaas_elec_prop(0.2, 0, 0, 300);
props3 = algaas_elec_prop(0.3, 0, 0, 300);
props4 = algaas_elec_prop(0.4, 0, 0, 300);

A0 = hbar^2/(2*props0.meff_e*9.11e-31)*(3*pi^2)^(2/3)*1e4/q; % the last two factors account for cm-3 to m-3 and J to eV. 
A1 = hbar^2/(2*props1.meff_e*9.11e-31)*(3*pi^2)^(2/3)*1e4/q;
A2 = hbar^2/(2*props2.meff_e*9.11e-31)*(3*pi^2)^(2/3)*1e4/q;
A3 = hbar^2/(2*props3.meff_e*9.11e-31)*(3*pi^2)^(2/3)*1e4/q;
A4 = hbar^2/(2*props4.meff_e*9.11e-31)*(3*pi^2)^(2/3)*1e4/q;

EF0 = A0*y.^(2/3);
EF1 = A1*y.^(2/3);
EF2 = A2*y.^(2/3);
EF3 = A3*y.^(2/3);
EF4 = A4*y.^(2/3);

figure;
% graphics = semilogx(y, EF0, y, EF1, y, EF2, y, EF3, y, EF4);
graphics = loglog(y, EF0, y, EF1, y, EF2, y, EF3, y, EF4);
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Carrier concentration cm^-^3','FontSize',11);
ylabel('Band filling induced shift E_F (eV)','FontSize',11);
legend('x = 0.0', 'x = 0.1', 'x = 0.2', 'x = 0.3', 'x = 0.4','Location','southeast');

%% Band tailing effect, equation (7) in the paper
B0 = (q^3*hbar/(4*3^(1/6)*pi^(1/3)*(props0.eps_r_static*eps0)^(3/2)*(9.11e-31*props0.meff_e)^0.5))^0.5*10^2.5/q;
B1 = (q^3*hbar/(4*3^(1/6)*pi^(1/3)*(props1.eps_r_static*eps0)^(3/2)*(9.11e-31*props0.meff_e)^0.5))^0.5*10^2.5/q;
B2 = (q^3*hbar/(4*3^(1/6)*pi^(1/3)*(props2.eps_r_static*eps0)^(3/2)*(9.11e-31*props0.meff_e)^0.5))^0.5*10^2.5/q;
B3 = (q^3*hbar/(4*3^(1/6)*pi^(1/3)*(props3.eps_r_static*eps0)^(3/2)*(9.11e-31*props0.meff_e)^0.5))^0.5*10^2.5/q;
B4 = (q^3*hbar/(4*3^(1/6)*pi^(1/3)*(props4.eps_r_static*eps0)^(3/2)*(9.11e-31*props0.meff_e)^0.5))^0.5*10^2.5/q;

ET0 = B0*y.^(5/12);
ET1 = B1*y.^(5/12);
ET2 = B2*y.^(5/12);
ET3 = B3*y.^(5/12);
ET4 = B4*y.^(5/12);

figure;
% graphics = semilogx(y, EF0, y, EF1, y, EF2, y, EF3, y, EF4);
graphics = loglog(y, ET0, y, ET1, y, ET2, y, ET3, y, ET4);
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Carrier concentration cm^-^3','FontSize',11);
ylabel('Band tailing induced shift E_T (eV)','FontSize',11);
legend('x = 0.0', 'x = 0.1', 'x = 0.2', 'x = 0.3', 'x = 0.4','Location','southeast');

%% Carrier-Carrier interaction effects
% in the above reference, no expression was found. IEEE J. Quantum
% Electronics Vol 26, Issue 1, Page 113, 1990 An xpression was found
% according equations (16), (17) and (18)

%Ref 2
EC0 = 0.125/props0.eps_r_static.*((y/(props0.chi_cr)-1).^(1/3)) .* heaviside(-(1-y/(props0.chi_cr)));
EC1 = 0.125/props1.eps_r_static.*((y/(props1.chi_cr)-1).^(1/3)) .* heaviside(-(1-y/(props1.chi_cr)));
EC2 = 0.125/props2.eps_r_static.*((y/(props2.chi_cr)-1).^(1/3)) .* heaviside(-(1-y/(props2.chi_cr)));
EC3 = 0.125/props3.eps_r_static.*((y/(props3.chi_cr)-1).^(1/3)) .* heaviside(-(1-y/(props3.chi_cr)));
EC4 = 0.125/props4.eps_r_static.*((y/(props4.chi_cr)-1).^(1/3)) .* heaviside(-(1-y/(props4.chi_cr)));

%Ref 1, modified to match Ref 2.
EC00 = 2.4e-8*y.^(1/3)*(props0.eps_r_static/props0.eps_r_static);
EC01 = 2.4e-8*y.^(1/3)*(props0.eps_r_static/props1.eps_r_static);
EC02 = 2.4e-8*y.^(1/3)*(props0.eps_r_static/props2.eps_r_static);
EC03 = 2.4e-8*y.^(1/3)*(props0.eps_r_static/props3.eps_r_static);
EC04 = 2.4e-8*y.^(1/3)*(props0.eps_r_static/props4.eps_r_static);

figure;
graphics = loglog(y, EC0, y, EC1, y, EC2, y, EC3, y, EC4, y, EC00, y, EC01, y, EC02, y, EC03, y, EC04);
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Carrier concentration cm^-^3','FontSize',11);
ylabel('Carrier-Carrier interaction E_C\prime (eV)','FontSize',11);
legend('x = 0.0', 'x = 0.1', 'x = 0.2', 'x = 0.3', 'x = 0.4','Location','southeast');

figure; % this is for checking Fig. 5 IEEE Journal of Quantum Electronics, Vol 26, No 1, P113, (1990)
graphics = semilogx(y, EC0/0.001);
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Carrier concentration (cm^-^3)','FontSize',13);
ylim([0 70]);
ylabel('-\DeltaE_g (meV)','FontSize',13);
legend('x = 0.0', 'x = 0.1', 'x = 0.2', 'x = 0.3', 'x = 0.4','Location','southeast');

%% Effective energy band gaps considering the above three effects

% Eg0 = props0.Eg + EF0 - ET0 -EC0;
% Eg1 = props1.Eg + EF1 - ET1 -EC1;
% Eg2 = props2.Eg + EF2 - ET2 -EC2;
% Eg3 = props3.Eg + EF3 - ET3 -EC3;
% Eg4 = props4.Eg + EF4 - ET4 -EC4;

Eg0 = props0.Eg + EF0 - ET0 -EC00;
Eg1 = props1.Eg + EF1 - ET1 -EC01;
Eg2 = props2.Eg + EF2 - ET2 -EC02;
Eg3 = props3.Eg + EF3 - ET3 -EC03;
Eg4 = props4.Eg + EF4 - ET4 -EC04;

figure;
graphics = semilogx(y, Eg0, y, Eg1, y, Eg2, y, Eg3, y, Eg4);
%graphics = loglog(y, Eg0, y, Eg1, y, Eg2, y, Eg3, y, Eg4);
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Carrier concentration cm^-^3','FontSize',11);
ylabel('Effective bandgaps E_g(eV)','FontSize',11);
legend('x=0', '0.1', '0.2', '0.3', '0.4');


%% check figure 2 of Ref [1]
El = 0.0077+3e-21*y;

figure;

graphics = semilogx(y,El);
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Carrier concentration cm^-^3','FontSize',11);
ylabel('E_l (eV)','FontSize',11);
xlim([min(y) max(y)]);
% ylim([10^0 10^5]);
% legend('N = 1E16', 'N = 5E16', 'N = 1E17', 'N = 5E17', 'N = 1E18','N = 5E18','Location','southeast');

%% Calculate the absorption coefficient

z = 0.2; % look at Aluminum concentration z
T = 300;
propsz = algaas_elec_prop(z, 0, 0, T); 


E = 0.8:0.001:2.0;
y = cat(2, linspace(1e17,5e17,2),linspace(1e18,5e18,2));

Az = hbar^2/(2*propsz.meff_e*9.11e-31)*(3*pi^2)^(2/3)*1e4/q;
Bz = (q^3*hbar/(4*3^(1/6)*pi^(1/3)*(propsz.eps_r_static*eps0)^(3/2)*(9.11e-31*propsz.meff_e)^0.5))^0.5*10^2.5/q;

EFz = NaN(1,length(y));
ETz = NaN(1,length(y));
ECz = NaN(1,length(y));
Egz = NaN(1,length(y));

mreduced = NaN(1,length(y));

El = NaN(1,length(y)); % equation (2) in Ref [1]
absorp = NaN(length(y),length(E));
absorp_original= NaN(length(y),length(E));

Es = 0.005; % make plot continuous at the boundary

for m = 1:length(y)
    
    mreduced(m) = 7.39e5*(2*((1-6.34e-14*y(m)^(2/3))/propsz.meff_e+1/propsz.meff_h)^(-1))^(3/2);
    
    El(m) = 0.0077+3e-21*y(m);
    EFz(m) = Az*y(m).^(2/3);    
    ETz(m) = Bz*y(m).^(5/12);
%     ECz(m) = 2.4e-8*y(m).^(1/3)*(props0.eps_r_static/propsz.eps_r_static);
    ECz(m) = 0.125/propsz.eps_r_static.*((y(m)/(propsz.chi_cr)-1).^(1/3)) .* heaviside(-(1-y(m)/(propsz.chi_cr)));
    Egz(m) = propsz.Eg + EFz(m)-ETz(m)-ECz(m);
    for n = 1:length(E)
        if (E(n)<=Egz(m))
            absorp(m,n) = 2000*exp((E(n)-Egz(m))/El(m));
        else
%             absorp(m,n) = 30000*(1-5.59e-14*y(m)^(2/3))^(-3/2)*(E(n)-Egz(m)+Es)^(1/2);
            absorp(m,n) = mreduced(m)*(E(n)-Egz(m)+Es)^(1/2);

        end
    end
    
        mreduced(m) = 7.39e5*(2*((1-6.34e-14*y(m)^(2/3))/propsz.meff_e+1/propsz.meff_h)^(-1))^(3/2);
    
    for n = 1:length(E)
        if (E(n)<=propsz.Eg)
            absorp_original(m,n) = 2000*exp((E(n)-propsz.Eg)/(kb*T/q));
%             absorp_original(m,n) = 0;
        else
%             absorp(m,n) = 30000*(1-5.59e-14*y(m)^(2/3))^(-3/2)*(E(n)-Egz(m))^(1/2);
            absorp_original(m,n) = (7.39e5*(2*(1/propsz.meff_e+1/propsz.meff_h)^(-1))^(3/2))*(E(n)-propsz.Eg+Es)^(1/2);
        end
    end
    
end

figure;
graphics = semilogy(E, absorp(1,:), E, absorp(2,:), E, absorp(3,:), E, absorp(4,:));
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Photon energy (ev)','FontSize',11);
ylabel('absorption coefficient (cm^-^1)','FontSize',11);
xlim([min(E) max(E)]);
ylim([10^-1 10^5]);
legend('N = 1E17', 'N = 5E17', 'N = 1E18','N = 5E18','Location','northwest');

figure;
graphics = semilogy(E, absorp_original(1,:), E, absorp_original(2,:), E, absorp_original(3,:), E, absorp_original(4,:));
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Photon energy (ev)','FontSize',11);
ylabel('original absorption coefficient (cm^-^1)','FontSize',11);
xlim([min(E) max(E)]);
% ylim([10^-1 10^5]);
legend('N = 1E17', 'N = 5E17', 'N = 1E18','N = 5E18','Location','northwest');

%% Calculate and plot defference
difference = absorp - absorp_original;

figure;
graphics = semilogy(E, -difference(1,:), E, -difference(2,:), E, -difference(3,:), E, -difference(4,:));
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Photon energy (ev)','FontSize',11);
ylabel('differential absorption coefficient (cm^-^1)','FontSize',11);
xlim([min(E) max(E)]);
% xlim([1.25 1.5]);
ylim([10^-1 10^5]);
legend('N = 1E17', 'N = 5E17', 'N = 1E18','N = 5E18','Location','northwest');

%% Kramers-Kronig relation

delta_n = NaN(length(y),length(E));

for k = 1:length(y)
    delta_n(k,:) = kkintegral(E,difference(k,:),E);
end

figure;
graphics = plot(E,-delta_n(1,:),E,-delta_n(2,:),E,-delta_n(3,:),E,-delta_n(4,:));
set(gca,'FontSize',11);
set(gca,'LineWidth',2);
set(graphics,'LineWidth',2);
xlabel('Photon energy (ev)','FontSize',11);
ylabel('\Deltan','FontSize',11);
xlim([min(E) max(E)]);
% ylim([10^1 10^5]);
legend('N = 1E17', 'N = 5E17', 'N = 1E18','N = 5E18','Location','northwest');