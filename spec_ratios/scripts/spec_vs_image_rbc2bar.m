%Aryil Bechtel 2021
%{
    plot rbc/bar ratio from spectroscopy vs rbc/bar in gas exchange
    reports

    import vars:
        -spec_bar2gas: bar/gas from spectroscopy
        -spec_rbc2gas: rbc/gas from spectroscopy
        -image_rbc2bar: rbc/bar from imaging
%}

clear;close all;

%load data
load('Duke_18_22/Xenon_lab/project_spec_ratios/data/19_Q4_image_rbc2bar.mat');
load('Duke_18_22/Xenon_lab/project_spec_ratios/data/19_Q4.mat');

%calc rbc2bar from spectroscopy ratios
spec_rbc2bar = spec_rbc2gas./spec_bar2gas;

%plot
figure();
x = linspace(0,max(max(spec_rbc2bar),max(image_rbc2bar))+0.1);
plot(x,x);
legend({'(rbc/bar)_{spec} = (rbc/bar)_{image}'},'AutoUpdate','off','EdgeColor','none')
hold on;
plot(spec_rbc2bar,image_rbc2bar,'*');
xlabel('(rbc/bar)_{spec}');
ylabel('(rbc/bar)_{image}');


