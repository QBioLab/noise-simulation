clear all
%         1   2     3     4      5    6   7    8    9      10    11  12     13     14     15     16   17   18     19     20     21     22     23     24   25     26
%param_=[P300 HDAC  Gt    k10    I0   n6  k_1  k20  k_20   A0    n2  k30    e31    k_30   k4     k_4  k5   k_5    k60    U0     e6     k7     k_7    k8   k9     k_9];
param  =[1e4  1,1e5 5.4e5 8.3e10 89   5   0.36 2e-3 2.5e-2 0.047 14  1.67e2 2.8e-3 4.5e-5 3.6e-7 4e-2 9e-5 2.2e-4 2.8e-2 1.1e-1 2.3e-4 1.1e-2 4.5e-2 1.67 1.3e-3 3.3e-4]
file_simu='AMV40-simu-result.mat';
[dt]=V40_AM_para_ic_rev3(file_simu,param)
