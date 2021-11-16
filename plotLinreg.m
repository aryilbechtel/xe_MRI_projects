function [mdl] = plotLinreg(in1,in2,xTitle,yTitle)

%{
plot scatter plot with linear regression model

in1 and in2 are arrays of same length
in1: indep var
in2: dep var 
xTitle: in1 label
yTitle: in2 label

inTitle: string title of graph
leg1 and leg2:  legend labels of points (leg1) and LSRL (leg2)
rSq is R^2 value of fit

%}

inTitle = '';
leg1 = '2019 Q4';
leg2 = 'linear fit';

%extract fit parameters
mdl = fitlm(in1,in2);
intercept = table2array(mdl.Coefficients(1,1));
slope = table2array(mdl.Coefficients(2,1));
rSq = mdl.Rsquared.Ordinary;

fprintf('r^{2}: ');
fprintf('%4.4f',rSq);
fprintf(newline);

fprintf('slope: ');
fprintf('%4.4f',slope);
fprintf(newline);

fprintf('intercept: ');
fprintf('%4.4f',intercept);

%plot
x1 = linspace(min(in1),max(in1),length(in1));
y1 = slope*x1 + intercept;

figure();
h1 = plot(in1,in2,'*');
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
hold on;
plot(x1,y1);
title(inTitle);
xlabel(xTitle);
ylabel(yTitle);
leg2 = append(leg2,' r^{2}: ',num2str(rSq));
legend(leg2);


figure();
plotResiduals(mdl,'fitted');
end

