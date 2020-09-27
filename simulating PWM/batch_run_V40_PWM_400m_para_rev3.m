clear all
%         1   2     3     4      5    6   7    8    9      10    11  12     13     14     15     16   17   18     19     20     21     22     23     24   25     26
%param_=[P300 HDAC  Gt    k10    I0   n6  k_1  k20  k_20   A0    n2  k30    e31    k_30   k4     k_4  k5   k_5    k60    U0     e6     k7     k_7    k8   k9     k_9];
param  =[1e4  1,1e5 5.4e5 8.3e10 89   5   0.36 2e-3 2.5e-2 0.047 14  1.67e2 2.8e-3 4.5e-5 3.6e-7 4e-2 9e-5 2.2e-4 2.8e-2 1.1e-1 2.3e-4 1.1e-2 4.5e-2 1.67 1.3e-3 3.3e-4]
period=400;
dt=1; tm=2800;  t=0:dt:tm; nT=length(t);
dtopen=[1 2 3 4 7 10 16 25 40 60 100 160 250 400]';
nM=length(dtopen); nPW=length(tm/400);
Ton=repmat([0:400:2400],nM,1); 
Toff=Ton+repmat(dtopen,1,nPW);
Ton=Ton+400;
nic=21; % number of initial values for chromosome accessibility and histone AC
I0=100; % light intensity for "ON" period of PWM, needs to be adjusted.
I1=0; % light intensity for "Off" period of PWM
I=repmat([I0 I1],1,tm/400);
file_simu='PWM400-simu-result.mat';
[dt]=V40_PWM_para_ic_rev3(file_simu,param,I,Ton, Toff,tm,dt,nic,l1(i1));
