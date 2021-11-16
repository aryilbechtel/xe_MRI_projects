%Aryil Bechtel 2021

%{
    Graph Bland-Altman plots for rbc/gas and bar/gas from spectroscopy and
    gas exchange imaging, color-coded for additional variables 
    (depending on section of code executed)

    import vars:
        -spec_bar2gas: bar/gas from spectroscopy
        -spec_rbc2gas: rbc/gas from spectroscopy
        -GX_bar_mean: mean bar/gas from imaging
        -GX_rbc_mean: mean rbc/gas from imaging
%}

%% 
%{
    additional import vars:
        -spec_linewidth: linewidth of dedicated gas spectral peak
%}

clear;close all;

colbarToggle = 0; %set colbarToggle = 1 for colorbar based on spec_linewidth
                  %set colbarToggle = 0 for no colorbar (normal
                  %Bland-Altman)

%load data
[file, path] = uigetfile('*.*', 'Select file');  % starts in current dir 
file_with_path = strcat(path, file);  % join path and filename to open
load(file_with_path);

%{
%remove data with issues
spec_bar2gas(totIssues) = [];
spec_rbc2gas(totIssues) = [];
GX_bar_mean(totIssues) = [];
GX_rbc_mean(totIssues) = [];
spec_linewidth(totIssues) = [];
%}

%declare arrays
varLength = length(spec_bar2gas);
rbcDiff = zeros(1,varLength);
rbcMean = zeros(1,varLength);
barDiff = zeros(1,varLength);
barMean = zeros(1,varLength);

for i=1:varLength
    
    %calc diff between ratios
    rbcDiff(i) = spec_rbc2gas(i) - GX_rbc_mean(i);
    barDiff(i) = spec_bar2gas(i) - GX_bar_mean(i);
    
    %calc mean of ratios
    rbcMean(i) = (spec_rbc2gas(i) + GX_rbc_mean(i))/2;
    barMean(i) = (spec_bar2gas(i) + GX_bar_mean(i))/2;
    
end

rbcDiff = (rbcDiff./rbcMean).*100;
barDiff = (barDiff./barMean).*100;

%mean and std of rbcDiff and barDiff
meanRbcDiff = mean(rbcDiff);
stdRbcDiff = std(rbcDiff);

meanBarDiff = mean(barDiff);
stdBarDiff = std(barDiff);


%rbc plot
figure();
if colbarToggle==1
    colormap cool;
    rbcC = colorbar;
    rbcC.Label.String = 'dedicated gas linewidth (Hz)';
    rbcC.Label.FontSize = 17;
    h1 = scatter(rbcMean,rbcDiff,30,spec_linewidth,'filled');
else
    h1 = scatter(rbcMean,rbcDiff,30,'filled');
end
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
xlabel('mean rbc/gas','FontSize',17);
ylabel('(rbc/gas)_{spec} - (rbc/gas)_{image}','FontSize',17);
hold on;
%hLineRbc = yline(0,'k');
%set(get(get( hLineRbc, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

xl = xlim;

yline(meanRbcDiff,'b');
%text(xl(1),meanRbcDiff,'mean','FontSize',12,'Color','r');

yline(meanRbcDiff+2*stdRbcDiff,'--');
%text(xl(1),(meanRbcDiff+2*stdRbcDiff),'mean + 2 SD','FontSize',12,'Color','r');

yline(meanRbcDiff-2*stdRbcDiff,'--');
%text(xl(1),(meanRbcDiff-2*stdRbcDiff),'mean - 2 SD','FontSize',12,'Color','r');

legend('mean','$\pm$ 2 SD','interpreter','latex','FontSize',17);
hold off;


%barrier plot
figure();
if colbarToggle==1
    colormap cool;
    barC = colorbar;
    barC.Label.String = 'dedicated gas linewidth (Hz)';
    barC.Label.FontSize = 17;
    h3 = scatter(barMean,barDiff,30,spec_linewidth,'filled');
else
    h3 = scatter(barMean,barDiff,30,'filled');
end
set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
xlabel('mean bar/gas','FontSize',17);
ylabel('(bar/gas)_{spec} - (bar/gas)_{image}','FontSize',17);
hold on;
%hLineBar = yline(0,'k');
%set(get(get( hLineBar, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

xl = xlim;

yline(meanBarDiff,'b');
%text(xl(1),meanBarDiff,'mean','FontSize',12,'Color','r');

yline(meanBarDiff+2*stdBarDiff,'--');
%text(xl(1),(meanBarDiff+2*stdBarDiff),'mean + 2 SD','FontSize',12,'Color','r');

yline(meanBarDiff-2*stdBarDiff,'--');
%text(0.5*xl(2),(meanBarDiff-2*stdBarDiff),'mean - 2 SD','FontSize',12,'Color','r');

legend('mean','$\pm$ 2 SD','interpreter','latex','FontSize',17);
hold off;

%histogram of rbc differences
numbins=6;
figure();
histogram(rbcDiff,numbins);
xlabel('(rbc/gas)_{spec} - (rbc/gas)_{image}','FontSize',14);
ylabel('counts','FontSize',14);

%histogram of barrier differences
figure();
histogram(barDiff,numbins);
xlabel('(bar/gas)_{spec} - (bar/gas)_{image}','FontSize',14);
ylabel('counts','FontSize',14);
%% 
%{
    additional import vars:
        -totIssues: logical array indicating which subjects had issues
        during scan (mostly from inhalation or breath hold)
%}

clear;close all;

%load data
load('Duke_18_22/Xenon_lab/project_spec_ratios/data/19_Q4.mat');
load('Duke_18_22/Xenon_lab/project_spec_ratios/data/19_Q4_issues.mat');

%declare arrays
varLength = length(spec_bar2gas);
rbcDiff = zeros(1,varLength);
rbcMean = zeros(1,varLength);
barDiff = zeros(1,varLength);
barMean = zeros(1,varLength);

%{
%try taking log of original data
spec_rbc2gas = log(spec_rbc2gas);
spec_bar2gas = log(spec_bar2gas);
GX_rbc_mean = log(GX_rbc_mean);
GX_bar_mean = log(GX_bar_mean);
%}

for i=1:varLength
    
    %calc diff between ratios
    rbcDiff(i) = spec_rbc2gas(i) - GX_rbc_mean(i);
    barDiff(i) = spec_bar2gas(i) - GX_bar_mean(i);
    
    %calc mean of ratios
    rbcMean(i) = (spec_rbc2gas(i) + GX_rbc_mean(i))/2;
    barMean(i) = (spec_bar2gas(i) + GX_bar_mean(i))/2;
    
end

%mean and std of rbcDiff and barDiff
meanRbcDiff = mean(rbcDiff);
stdRbcDiff = std(rbcDiff);

meanBarDiff = mean(barDiff);
stdBarDiff = std(barDiff);

%rbc plot
figure();
colormap cool;
scatter(rbcMean,rbcDiff,30,totIssues,'filled');
xlabel('mean rbc/gas','FontSize',14);
ylabel('(rbc/gas)_{spec} - (rbc/gas)_{image}','FontSize',14);
hold on;
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'. magenta','MarkerSize',30);
h(2) = plot(NaN,NaN,'. cyan','MarkerSize',30);
legend(h, 'issues','no issues','AutoUpdate','off');
yline(0,'k');

yline(meanRbcDiff,'b');
text(min(rbcMean),1.2*meanRbcDiff,'mean','FontSize',14);

yline(meanRbcDiff+2*stdRbcDiff,'--');
text(min(rbcMean),1.05*(meanRbcDiff+2*stdRbcDiff),'mean + 2 SD','FontSize',14);

yline(meanRbcDiff-2*stdRbcDiff,'--');
text(min(rbcMean),0.90*(meanRbcDiff-2*stdRbcDiff),'mean - 2 SD','FontSize',14);
hold off;


%barrier plot
figure();
colormap cool;
scatter(barMean,barDiff,30,totIssues,'filled');
xlabel('mean bar/gas','FontSize',14);
ylabel('(bar/gas)_{spec} - (bar/gas)_{image}','FontSize',14);
hold on;
h(1) = plot(NaN,NaN,'. magenta','MarkerSize',30);
h(2) = plot(NaN,NaN,'. cyan','MarkerSize',30);
legend(h, 'issues','no issues','AutoUpdate','off');
yline(0,'k');

yline(meanBarDiff,'b');
text(min(barMean),1.3*meanBarDiff,'mean','FontSize',14);

yline(meanBarDiff+2*stdBarDiff,'--');
text(min(barMean),1.08*(meanBarDiff+2*stdBarDiff),'mean + 2 SD','FontSize',14);

yline(meanBarDiff-2*stdBarDiff,'--');
text(min(barMean),0.87*(meanBarDiff-2*stdBarDiff),'mean - 2 SD','FontSize',14);

hold off;

%histogram of rbc differences
numbins=6;
figure();
histogram(rbcDiff,numbins);
xlabel('(rbc/gas)_{spec} - (rbc/gas)_{image}','FontSize',14);
ylabel('counts','FontSize',14);

%histogram of barrier differences
figure();
histogram(barDiff,numbins);
xlabel('(bar/gas)_{spec} - (bar/gas)_{image}','FontSize',14);
ylabel('counts','FontSize',14);

%%
%{
    import vars:
        -barDiff: difference between spec bar/gas and GX bar/gas
        -rbcDiff: difference between spec rbc/gas and GX rbc/gas
        -barMean: mean of spec bar/gas and GX bar/gas
        -rbcMean: mean of spec rbc/gas and GX rbc/gas
        -spec_linewidth: linewidth of dedicated gas spectral peak
%}

clear;close all;

%load data
load('Duke_18_22/Xenon_lab/project_spec_ratios/data/19_Q4_filtered_add_subs/19_Q4_filtered_add_subs_BlandAlt.mat');

%mean and std of rbcDiff and barDiff
meanRbcDiff = mean(rbcDiff);
stdRbcDiff = std(rbcDiff);

meanBarDiff = mean(barDiff);
stdBarDiff = std(barDiff);

%rbc plot
figure();
colormap cool;
scatter(rbcMean,rbcDiff,30,spec_linewidth,'filled');
rbcC = colorbar;
rbcC.Label.String = 'dedicated gas linewidth (Hz)';
rbcC.Label.FontSize = 17;
xlabel('mean rbc/gas','FontSize',17);
ylabel('(rbc/gas)_{spec} - (rbc/gas)_{image}','FontSize',17);
hold on;
yline(0,'k');

xl = xlim;

yline(meanRbcDiff,'b');
%text(xl(1),meanRbcDiff,'mean','FontSize',12,'Color','r');

yline(meanRbcDiff+2*stdRbcDiff,'--');
%text(xl(1),(meanRbcDiff+2*stdRbcDiff),'mean + 2 SD','FontSize',12,'Color','r');

yline(meanRbcDiff-2*stdRbcDiff,'--');
%text(xl(1),(meanRbcDiff-2*stdRbcDiff),'mean - 2 SD','FontSize',12,'Color','r');
hold off;

%{
%rbc percentage plot
figure();
colormap cool;
scatter(rbcMean,(rbcDiff./rbcMean).*100,30,spec_linewidth,'filled');
rbcC = colorbar;
rbcC.Label.String = 'dedicated gas linewidth (Hz)';
xlabel('mean rbc/gas');
ylabel('diff/mean (%)');
hold on;
yline(0,'k');
hold off;
%}

%barrier plot
figure();
colormap cool;
scatter(barMean,barDiff,30,spec_linewidth,'filled');
barC = colorbar;
barC.Label.String = 'dedicated gas linewidth (Hz)';
barC.Label.FontSize = 17;
xlabel('mean bar/gas','FontSize',17);
ylabel('(bar/gas)_{spec} - (bar/gas)_{image}','FontSize',17);
hold on;
yline(0,'k');

xl = xlim;

yline(meanBarDiff,'b');
%text(xl(1),meanBarDiff,'mean','FontSize',12,'Color','r');

yline(meanBarDiff+2*stdBarDiff,'--');
%text(xl(1),(meanBarDiff+2*stdBarDiff),'mean + 2 SD','FontSize',12,'Color','r');

yline(meanBarDiff-2*stdBarDiff,'--');
%text(0.5*xl(2),(meanBarDiff-2*stdBarDiff),'mean - 2 SD','FontSize',12,'Color','r');

hold off;

%{
%barrier percentage diff plot
figure();
colormap cool;
scatter(barMean,(barDiff./barMean).*100,30,spec_linewidth,'filled');
barC = colorbar;
barC.Label.String = 'dedicated gas linewidth (Hz)';
xlabel('mean bar/gas');
ylabel('diff/mean (%)');
hold on;
yline(0,'k');
hold off;
%}

%histogram of rbc differences
numbins=6;
figure();
histogram(rbcDiff,numbins);
xlabel('(rbc/gas)_{spec} - (rbc/gas)_{image}','FontSize',17);
ylabel('counts','FontSize',17);

%histogram of barrier differences
figure();
histogram(barDiff,numbins);
xlabel('(bar/gas)_{spec} - (bar/gas)_{image}','FontSize',17);
ylabel('counts','FontSize',17);

