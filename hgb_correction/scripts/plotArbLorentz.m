%{
Script to plot demo lorentzian distributions for gas, RBC, and barrier
spectra

Aryil Bechtel 2021
%}

clear; close all;

ppm = [-50:0.1:250];
gasPos = 0;
barPos = 197;
rbcPos = 218;
hwhmGas = 2;
hwhmBar = 10;
hwhmRbc = 30;
hGas = 1/2;
hBar = 1;
hRbc = 1/2;
hRbc2 = 1/5;

gasLorentz = arbitraryLorentz(ppm, gasPos, hwhmGas, hGas);
barLorentz = arbitraryLorentz(ppm, barPos, hwhmBar, hBar);
rbcLorentz = arbitraryLorentz(ppm, rbcPos, hwhmRbc, hRbc);
rbcLorentz2 = arbitraryLorentz(ppm, rbcPos, hwhmRbc, hRbc2);
figure();
plot(ppm, gasLorentz, 'LineWidth', 3, 'Color', '#0072BD');
hold on;
plot(ppm, barLorentz, 'LineWidth', 3,'Color','#77AC30');
plot(ppm, rbcLorentz, 'LineWidth', 3, 'Color', '#D95319');
plot(ppm, rbcLorentz2, 'LineWidth', 3, 'Color', '#D95319', 'LineStyle','--');
set(gca,'FontSize',18);
set(gca, 'XDir','reverse')
ylim([0,1.1])
set(gca,'YTickLabel',[]);

legend('gas','barrier','high RBC','low RBC');
xlabel('chemical shift (ppm)')
ylabel('absolute magnitude')
