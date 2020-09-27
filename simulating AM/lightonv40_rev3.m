function dydt = lightonv40_rev3(t,y,flag,para) %
[n6,k1,Gt,k_1, k20,k_20,A0,n2, k30,e31, k_30,HDC, k40,k_40, k50,P300, k_50, k60,U0, e6,k7, k_7,k8, k9, k_9]=para{:};
dydt=zeros(9,1);

% y(1) : [GAVPO dimer] [G2], free GAVPO [G]=[Gt]-2G2;
% k1=k1_0*I^n/(I^1.5+I0^1.5) calculated from outside, contains light
dydt(1) = k1* (Gt-2*y(1))* (Gt-2*y(1)) - k_1*y(1); 
% y(2) : CO chromatin opened states of two nuecleasomes CC=2-CO
    %HAC=y(3);
dydt(2) = k20* (2-y(2)) - k_20*y(2)*A0^n2/(A0^n2+y(3)^n2); 
% y(3) : HAC Histone acetylation states H=2-HAC 
    %UGB=y(5); 
dydt(3) = k30*(y(5)+e31)* (2-y(3)) - k_30*HDC*y(3);
% y(4) : UAS bounded with GAVPO2: UU=[2.5*CO]-UB
    %CO25=(2.5*y(2); 
dydt(4) = k40 * y(1) *  (2.5*y(2)-y(4)) - k_40*y(4);
% y(5) : UAS-G2 bounded with p300: UGU = UB - UGB, 
   % UB=y(4); 
dydt(5) = k50*P300* (y(4)-y(5)) - k_50*y(5);
% y(6) : nuclear mRNA transcription from P300 occupied promoter (Pt=8), with synergism (hill coef n1) decrease due to translocation out of nuclear
    %UGB=y(5);%floor(y(5)+0.5); % descretize UAS available states from ODE
dydt(6) = 8*k60*(y(5)^n6/(y(5)^n6+U0^n6)+e6) - k7*y(6);
% y(7):  cytosolic mRNA, increase due to translocation, decrease due to degradation
dydt(7) = k7*y(6) - k_7*y(7);
% y(8): immature FP: increase by translation, decrease by maturation
dydt(8) = k8*y(7) - k9*y(8);
% y(9): mature FP , increase by maturation, decrease by degradation
dydt(9) = k9*y(8) - k_9*y(9);