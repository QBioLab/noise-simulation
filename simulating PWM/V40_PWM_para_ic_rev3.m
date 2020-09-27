function [dt]=V40_PWM_para_ic_rev3(file_simu,param,I,Ton, Toff,tm,dt,nic,simu_no) %i1: the first simulation, with no parameter permutation
    npara=26;
    np=15;
    % permutate HDAC Gt, 2^(-7:7) with other set to const
    perG=2.^([-7:7]'); % to minic GAVPO variations
    perH=2.^(-7:7);% to minic HDAC variations or LMK inhibition
    perG2=repmat(perG,1,np);
    perH2=repmat(perH,np,1);
    per=[perH2(:) perG2(:)];
    para3=[ones(np*np,1) per]; 
    npa3=length(para3(:,1));
    param0=[para3 ones(npa3,23)];
    param1=repmat(param,npa3,1).*param0;

    t=0:dt:tm; nT=length(t); 
    nv=9;
    [nM nPW]=size(Ton);
    % save(file_simu,'param','param0','param1');

    Tpwm=ones(nM,tm);
    for i1=1:nM
        for i2=1:nPW
            Tpwm(i1,Toff(i1,i2)+1:Ton(i1,i2))=0;
        end
    end
    I0=Tpwm*I(1);
    % set initial chromosome accessibility and histone AC levels
    [ic]=initial_conditions(nic,9);
    realn=zeros(nT,np*np,nM*nic);

    paramtotal=param1;
    out34=zeros(nT,nv,np*np,nM*nic);
    parfor i12=1:nic*nM
        [out12,realn12]=V40_para_per_10min_rev(i12,nM,ic,nT,nv,tm,I0,paramtotal);
        out34(:,:,:,i12)=out12;% time, variable, PWM, ic, 
        realn(:,:,i12)=realn12;
    end
     result.I=I0;
     result.t=t;
     result.param0=param;
     result.param1=param1;
     result.realnum=realn;
     result.ic=ic;
     result.outsIC=out34;
     save(file_simu,'param','param0','param1','result','-v7.3');
         tt=toc;    
         [simu_no tt/1000]    
end

function [ic]=initial_conditions(nic,nv)
    ic=zeros(nic,nv);% initial conditions
    % initialize G2
    ic(1:nic,1)=0;                
    % initialize nucleosome and histone-AC
    for i4=2:3
        dic=2/(nic-1);
        ic(1:nic,i4)=0+(0:dic:2);
    end
    ic(1:nic,4:nv)=0;  % light,initial C, variables, 100       
end




