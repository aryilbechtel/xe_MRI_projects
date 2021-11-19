%{
Script to create simple fits to MOXE hemoglobin correction functions

Aryil Bechtel 2021
%}

%% Generate correction functions
%clear all; close all;

%set MOXE constants
flip_angle = 20; % flip angle degrees for equivalent RBC/barrier analysis
hgb_ref = 15.0; % hgb g/dl for equivalent RBC/barrier analysis
hgb_deviation = 10;
TR = 15; % repetition time (ms) for equivalent RBC/barrier analysis

D = 3.3e-10; % tissue diffusion coefficient in m2/sec, old Ruppert measurement
lambda = 0.2; % general solubility of xenon in stuff - REVISIT this. Has implications for estimated RBC/gas and barrier/gas
d = 1e-5; % total septal wall thickness in meters (9.2 microns)
delta = 0.6e-6;%1.0e-6; % barrier thickness in m (literature is 0.6 microns, but use 1 micron)
t_x = 1.6; % capillary transit time in seconds. Chang uses 1.3 (a bit long, unless this is from right to left heart)
SA_to_V = 25000; % surface to volume in m^-1 (210+-50 cm^-1)
b = lambda*SA_to_V*d/2;  % signal normalization constant
T_exch=d^2/(pi^2*D); % principal exchange time constant (s) per eq 5 of MOXE. Chang quotes 26±13 in MRM 2014
nterms = 5; % number of terms in series summations; MOXE paper uses 5 terms

%evaluate correction functions
TR_equiv_eval = (TR/1000)/(1-cos(pi*flip_angle/180));

hgb_array_length = 60;
hgb_array = linspace(hgb_ref - hgb_deviation,hgb_ref + hgb_deviation,hgb_array_length); % array of hgb values to evaluate
rbc2barFactor_array = zeros(1,length(hgb_array));
rbcFactor_array = zeros(1,length(hgb_array));
barFactor_array = zeros(1,length(hgb_array));

for i = 1:length(rbc2barFactor_array)
    [rbc2barFactor_array(i), rbcFactor_array(i), barFactor_array(i)] = ...
        correctHgbMOXE(b,delta,d,t_x,T_exch,nterms,TR_equiv_eval,hgb_ref,hgb_array(i));
end

%% create fits
%polynomial fits
n = 3; %degree of polynomial
[r2bPoly3,r2bPoly3_errors] = polyfit(hgb_array,rbc2barFactor_array,n);
[rPoly3,rPoly3_errors] = polyfit(hgb_array,rbcFactor_array,n);
[bPoly3,bPoly3_errors] = polyfit(hgb_array,barFactor_array,n);

%exponential fits
expFit = fit(hgb_array',rbc2barFactor_array','exp1');

%% plot
close all;
figure();
plot(hgb_array,rbc2barFactor_array,'Color','k','LineWidth',3);
hold on;
plot(hgb_array,rbcFactor_array,'Color','#D95319','LineWidth',3);
plot(hgb_array,barFactor_array,'Color', '#77AC30','LineWidth',3);
%plot(hgb_array,polyval(r2bPoly3,hgb_array),'Color','r','LineStyle','--','LineWidth',4);
%plot(hgb_array,polyval(rPoly3,hgb_array),'Color','r','LineStyle','--','LineWidth',4);
%plot(hgb_array,polyval(bPoly3,hgb_array),'Color','r','LineStyle','--','LineWidth',4);
yline(1,'k');

a = sprintf(['Bar, RBC, and RBC/bar relative to Hgb=%2.1f g/dl (flip= %2.0f' char(176)...
    ', TR=%2.0f ms)'],hgb_ref,flip_angle,TR);
title(a);
xlabel('Hgb (g/dL)');
ylabel(texlabel('scaling factor (beta_Hgb)'));
xlim([10,20])

r2bEq1 = num2str(r2bPoly3(1),'%2.2e');
r2bEq2 = num2str(r2bPoly3(2),'%2.2e');
r2bEq3 = num2str(r2bPoly3(3),'%2.2e');
r2bEq4 = num2str(r2bPoly3(4),'%2.2e');
r2bEqtot = strcat(texlabel('beta_{Hgb,RBC/bar} = ('),r2bEq1,')Hgb^3 + (',r2bEq2,')Hgb^2 + (',...
    r2bEq3,')Hgb + (',r2bEq4,')');

rEq1 = num2str(rPoly3(1),'%2.2e');
rEq2 = num2str(rPoly3(2),'%2.2e');
rEq3 = num2str(rPoly3(3),'%2.2e');
rEq4 = num2str(rPoly3(4),'%2.2e');
rEqtot = strcat(texlabel('beta_{Hgb,RBC} = ('),rEq1,')Hgb^3 + (',rEq2,')Hgb^2 + (',...
    rEq3,')Hgb + (',rEq4,')');

bEq1 = num2str(bPoly3(1),'%2.2e');
bEq2 = num2str(bPoly3(2),'%2.2e');
bEq3 = num2str(bPoly3(3),'%2.2e');
bEq4 = num2str(bPoly3(4),'%2.2e');
bEqtot = strcat(texlabel('beta_{Hgb,bar} = ('),bEq1,')Hgb^3 + (',bEq2,')Hgb^2 + (',...
    bEq3,')Hgb + (',bEq4,')');

legend('RBC/Bar','RBC','bar','AutoUpdate','off');
text(10.3,1.72,r2bEqtot,'FontSize',14,'Color','k','FontWeight','bold');
text(10.3,1.62,rEqtot,'FontSize',14,'Color','#D95319','FontWeight','bold');
text(10.3,1.52,bEqtot,'FontSize',14,'Color','#77AC30','FontWeight','bold');