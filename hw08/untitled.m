close all;
xdata = [3 10 20 30 40 50];
iter1 = [17 35 67 107 137 50000];
iter2 = [20 249 909 1939 3262 4974];
figure 
semilogy(xdata,iter1,'-r*');
hold on
semilogy(xdata,iter2,'-b*');
legend('ShiftedQR','QR');
xlabel('Matrix Dimension','FontSize',12,'FontWeight','bold');
ylabel('Iterarions','FontSize',12,'FontWeight','bold');
title('Matrix Dimension vs Iterations','FontSize',12,'FontWeight','bold');


figure;
timeShiftedQR = [0 0.005 0.021 0.099 0.282 212.621];
timeQR = [0
0.012
0.256
1.793
7.379
22.997
 ];
timeQR = timeQR';
semilogy(xdata,timeShiftedQR,'-r*');
hold on;
semilogy(xdata,timeQR,'-b*');
legend('ShiftedQR','QR');
xlabel('Matrix Dimension','FontSize',12,'FontWeight','bold');
ylabel('CPU time(s)','FontSize',12,'FontWeight','bold');
title('Matrix Dimension vs CPU time','FontSize',12,'FontWeight','bold');
hold off


figure
lambda1 = [6.372281
67.840399
270.495189
608.253606
1081.115447
1689.080688
];
lambda2 = [2
20.431729
81.223819
182.544889
324.394506
506.772618
];
lambda3 = [0.627719
4.455992
17.235222
38.53868
68.364136
106.711318
];
lambdan = [0.627719
0.512543
0.503097
0.501373
0.500772
0.500494
];
lambdan_1 = [2
0.55164
0.512479
0.505511
0.503096
0.501978
];
lambdan_2 = [6.372281
0.629808
0.528819
0.512543
0.507004
0.500468
];
semilogy(xdata,lambda1,'-r*');
hold on 
semilogy(xdata,lambda2,'-g*');
semilogy(xdata,lambda3,'-b*');
semilogy(xdata,lambdan,'-c*');
semilogy(xdata,lambdan_1,'-m*');
semilogy(xdata,lambdan_2,'-k*');
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_n','\lambda_{n-1}','\lambda_{n-2}');
xlabel('Matrix Dimension','FontSize',12,'FontWeight','bold');
ylabel('Eigenvalue','FontSize',12,'FontWeight','bold');
title('Matrix Dimension vs Eigenvalues','FontSize',12,'FontWeight','bold');

