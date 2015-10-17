close all;
%timeLU=[0.026 0.156 1.219 9.895 78.276];
Vne = [0.482759 0.441080 0.407838 0.392824 0.382576 0.379957 0.373304    ];
Vsw = [0.275862 0.280752 0.294016 0.301635 0.307032 0.308423 0.311969    ];
Vse = [0.172414 0.211722 0.248441 0.265600 0.277288 0.280271 0.287845    ];

rperside = [2 4 10 20 40 50 100];
iter = [9 20 49 96 184 226 426];
Req = [1208.333333 850.525689 491.202658 307.389238 185.644983 156.853212 91.476893];
mat_size=[ 2 4 10 20];
time6 = [0 0.000 0.017 0.401 19.540 64.904 2005.768];
time4 = [0 0.002 0.012 0.510 28.149 104.032 6229.428];
%n_2=[78.276/4096,78.276/512,78.276/64,78.276/8,78.276];
figure
plot(rperside,iter);
xlabel('Resistors Per Side','FontSize',12,'FontWeight','bold');
 ylabel('Iterations','FontSize',12,'FontWeight','bold');
figure
%loglog(mat_size, timeLU,'-.r*');
 loglog(rperside, time6,'-.r*');
hold on
%timeCho=[0.031 0.218 0.877 5.612 54.396];
loglog(rperside, time4,':bs');

% figure
% plot(rperside, iter);
% loglog(mat_size, Vne,'-gs');
legend('HW06','HW04');
xlabel('Resistors Per Side','FontSize',12,'FontWeight','bold');
ylabel('Time(s)','FontSize',12,'FontWeight','bold');
%log(0.156/0.026)
grid on
hold off

figure;
plot(rperside, Vne);
hold on
plot(rperside, Vsw);
plot(rperside, Vse);

legend('Vne','Vsw','Vse');
xlabel('Resistors Per Side','FontSize',12,'FontWeight','bold');
ylabel('Voltage(V)','FontSize',12,'FontWeight','bold');
hold off
figure 
plot(rperside,Req);
xlabel('Resistors Per Side','FontSize',12,'FontWeight','bold');
ylabel('Req(Ohm)','FontSize',12,'FontWeight','bold');


 