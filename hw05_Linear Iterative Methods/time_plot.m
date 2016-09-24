time1=[0 0 0.026 0.156 1.219 9.895 78.276];
mat_size=[3 10 100 200 400 800 1600];
figure

xlabel('MATsize');
ylabel('Time(s)');
n_3=[nthroot(78.276,21),nthroot(78.276,18),nthroot(78.276,15),nthroot(78.276,12),nthroot(78.276,9),nthroot(78.276,3),78.276];

loglog([mat_size;mat_size], [time1;n_3], '-s')
lagend('LU Decomposition','n^3');
grid on
