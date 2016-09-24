close all;
xdata = [2^2     4^2    10^2    20^2    40^2];
ydata = [0.0010    0    0.0230    0.6360   31.5630];
figure 
loglog(xdata,ydata);
xlabel('Resistors Per Side','FontSize',12,'FontWeight','bold');
ylabel('CPU Time(s)','FontSize',12,'FontWeight','bold');
title('Analysis of CPU time vs resistor per side');
figure;

lambdaMin = [0.000541 0.000277 0.000094 0.000040 0.000017];
lambdaMax = [1 1 1 1 1];
loglog(xdata,lambdaMax);
hold on;
loglog(xdata,lambdaMin);
xlabel('Resistors Per Side','FontSize',12,'FontWeight','bold');
ylabel('Values of eigenvalues','FontSize',12,'FontWeight','bold');
legend('lambda_1','lambda_n');
title('Analysis of eigenvalue vs resistor per side');
hold off
figure
cond = [1849.470024 3607.461004 10653.368024 25212.310283 59700.888558];
loglog(xdata,cond);
xlabel('Resistors Per Side','FontSize',12,'FontWeight','bold');
ylabel('Values of condition numbers','FontSize',12,'FontWeight','bold');
title('Analysis of condition number vs resistor per side');


