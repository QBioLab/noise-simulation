function [out0,flag]=permlighton_V40_specific_params_rev3(param,I,ic,tm,dt)
% permparam: permutation of parameters, length: parameter no,
% % parameter=param_*permroot.^(permparam)
% t: unit min
% concentration: unit no./nuclear 1~3.3pM
% chromosome unit: local no.
% concentration of proteins controlling transcription, acetylation and deacetylation 
% param=[param3 param0];

% 26 parameters including 3 protein concentrations P300, HDAC, Gt, 
% param_=[P300_ HDAC_ Gt_ k10_ I0_ n_ k_1_ k20_ k_20_ A0_ n2_ k30_ e31_ k_30_  ...
%             k40_ k_40_ k50_ k_50_ k60_ U0_ e6_ k7_ k_7_ k8_ k9_  k_9_];
para=num2cell(param);
[P300, HDAC, Gt, k10, I0, n6, k_1, k20, k_20, A0, n2, k30, e31, k_30,  ...
            k40, k_40, k50, k_50, k60, U0, e6, k7, k_7, k8, k9,  k_9]=para{:};

% ic=zeros(1,9); 
options=[];
t=0:dt:tm; nt=length(t);
out0=zeros(nt,9); % time point, variables
k1=k10*I^1.5/(I^1.5+I0^1.5);
    parameter0=[n6 k1 Gt k_1  k20 k_20 A0 n2 k30 e31  k_30 HDAC  k40 k_40  k50 P300  k_50  k60 U0  e6 k7  k_7 k8  k9  k_9];
    para=num2cell(parameter0);
    [t y]=ode15s('lightonv40_rev3',0:dt:tm,ic,options,para);
     flag=1;
% the ODE model is a stiff problem. 
% Improper parameters could cause pathological solutions, 
% such as negative, out of bound or imaginary values
   if exist('y') 
        if length(y(:,1))==nt
            if isreal(y)
                % preventing overshooting or negative results in simulation
                y(find(y<0))=0; 
                y(find(y(:,4)>5),4)=5;
                y(find(y(:,5)>5),5)=5;
                out0=y;
            else
                % remove results with imaginary value
                out0(:,:)=1e-12;
                flag=0;
            end
        else
               % remove results with "blow up" scenario
               out0(:,:)=1e-12;
                flag=0;
        end
   else
       % remove results with "blow up" scenario
           out0(:,:)=1e-12;
           flag=0;
   end

