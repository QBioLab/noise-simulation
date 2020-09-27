function [dt]=V40_AM_para_ic_rev3(file_save,param_basic)
tic
w = warning ('off','all');
%        1    2    3  4   5   6  7   8   9    10 11 12  13  14   15  16   17  18   19  20 21 22 23  24 25  26      
%param_=[P300 HDAC Gt k10 I0  n6 k_1 k20 k_20 A0 n2 k30 e31 k_30 k40 k_40 k50 k_50 k60 U0 e6 k7 k_7 k8 k9  k_9];
npara=length(param_basic);
% time 
dt=100; tm=2800; t=0:dt:tm; nT=length(t);
% light intensity scanning
I=10.^([-1.4:0.2:2.4]); nI=length(I);
% permutate HDAC Gt, 2^(-7:7) with other const
perG=2.^([-7:7]'); % to minic GAVPO variations
perH=2.^(-7:7);% to minic HDAC variations or LMK inhibition
perG2=repmat(perG,1,15);
perH2=repmat(perH,15,1);
per=[perH2(:) perG2(:)];
para3=[ones(225,1) per]; 
npa3=length(para3(:,1));
param0=[para3 ones(npa3,23)];
% initial value finner
nic=21; 

ic=zeros(nI,nic,9,npa3);% initial conditions
outsIC=zeros(nT,9,nI,nic,npa3);
realnum=zeros(nic,nI,npa3);
param=param_basic.*param0;
for i2=1:npa3
    for i3=1:nI
%             ic(i3,1:nic,1,i2)=0;   
            for i4=2:3
                dic=2/(nic-1);
                ic(i3,1:nic,i4,i2)=0+(0:dic:2);
            end
            ic(i3,1:nic,4:9,i2)=0;  % light,initial C, variables, 100       
     end
end 
parfor i1=1:npa3       
     paramtotal=param(i1,:);
     ic1=ic(:,:,:,i1);
     [out34, realn]=V40_nI_iC(paramtotal,I,ic1,tm,dt,nT,nI,nic);
     outsIC(:,:,:,:,i1)=out34;
     realnum(:,:,i1)=realn;
end
     outsICs=permute(outsIC,[4 3 1 2 5]); % initial condition, light, time, variables, npa3
result.I=I
result.t=t;
result.param=param;
result.realnum=realnum;
result.ic=ic;
result.outsIC=outsICs;

%filename=sprintf('%sode15s2-v312-scan-ic-%04d-%04d.mat',folder,floor(filenum),i0);
save(file_save,'param_basic','param','result','-v7.3');
tt=toc;
temp=find(realnum(:)==0);     
[length(temp)/nI/nic tt/1000]  
end

function [out34, realn]=V40_nI_iC(paramtotal,I,ic1,tm,dt,nT,nI,nic)
    out34=zeros(nT,9,nI,nic);
    realn=zeros(nic,nI);
    flag1=1;
    i3=1;
    while i3<= nI && flag1 ==1
        i4=1;
        while i4<=nic && flag1 ==1
            ic0=ic1(i3,i4,:);
            ic0=ic0(:);            
            [out,flag]=permlighton_V40_specific_params(paramtotal,I(i3),ic0,tm,dt);
            out34(:,:,i3,i4)=out; % time, variable, Intensity, ic, para3
            realn(i4,i3)=flag;
            flag1=flag;
            i4=i4+1;
        end
        i3=i3+1;
    end
end
