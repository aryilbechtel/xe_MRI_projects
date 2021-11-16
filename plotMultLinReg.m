%perform multiple linear regression with model: 2 inputs, 1 output

y=true_refVolt_tot'; %response variable
x1=weight_tot'; %input var 1
x2=Href_amp_tot'; %input var 2
X=[ones(size(x1)) x1 x2];

%create fit coefficients
%b = regress(y,X);
[b,bint,r,rint,stats]=regress(y,X);

%plot fit
figure();
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):2:max(x1);
x2fit = min(x2):10:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
mesh(X1FIT,X2FIT,YFIT);
view(50,10);
hold off;

title('weight vs 1H ref vs true Xe ref');
xlabel('weight (kg)');
ylabel('1H ref (V)');
zlabel('true Xe ref (V)');

leg1='2019 Q2-4 and 2020 Q1,3-4';
leg2 = append('linear fit r^{2}: ',num2str(stats(1)));
legend(leg1,leg2);

%plot fit values vs residuals
fitvals = zeros(1,length(x1));
for i=1:length(x1)
    fitvals(i) = b(1) + b(2)*x1(i) + b(3)*x2(i);
end
figure();
plot(fitvals,r,'*');
hold on;
yline(0,'k');
title('fit vals v residuals');
xlabel('fit vals');
ylabel('residuals');


