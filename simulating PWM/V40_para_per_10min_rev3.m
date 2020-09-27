function [out12,realn12]=V40_para_per_10min_rev3(i12,nM,ic,nT,nv,tm,I0,paramtotal)
    i1=ceil(i12/nM);
    i2=i12-(i1-1)*nM;
    ic0=ic(i1,:);
    ic0=ic0(:);
    tm0=10;
    temp=reshape(I0(i2,:),10,280);
    npara3=length(paramtotal(:,1));
    out12=zeros(nT,nv,npara3);
    realn12=zeros(nT, npara3);
    for i1=1:npara3
        param=paramtotal(i1,:);
            for i3=10:10:tm % time pieces
                if i3==10                    
                    ic01=ic0;
                    ic01=ic01(:);
                    realn12(1,i1)=1;
                else
                    ic01=out12(i3-9,:,i1); % new IC from old end point
                    ic01=ic01(:);
                end
                out12(i3+1-10,:,i1)=ic01';
                temp=I0(i2,1:10)/I0(1,1);
                if sum(temp)==10 | sum(temp)==0
                    IPW=I0(i2,i3-9);   
                    tmm=10;
                    dt=0.1;
                    [out,flag]=permlighton_V40_specific_params_rev3(param,IPW,ic01,tmm,dt);
                    out12(i3-9+(1:10),:,i1)=out(11:10:end,:);% time, variable, PWM, ic, 
                    realn12(i3-9+(1:10),i1)=flag;
                else
                    tm1=sum(temp);
                    dt=0.1;
                    IPW=I0(i2,i3+1-10); 
                    [out,flag]=permlighton_V40_specific_params_rev3(param,IPW,ic01,tm1,dt);
                    out12(i3-9+(1:tm1),:,i1)=out(11:10:end,:);% time, variable, PWM, ic, 
                    realn12(i3-9+(1:tm1),i1)=flag;
                    
                    ic01=out(end,:); % new IC from old end point
                    ic01=ic01(:);
                    tm2=10-tm1;
                    IPW=I0(i2,i3+1-10+tm1); 
                    [out,flag]=permlighton_V40_specific_params_rev3(param,IPW,ic01,tm2,dt);
                    out12(i3-9+tm1+(1:tm2),:,i1)=out(11:10:end,:);% time, variable, PWM, ic, 
                    realn12(i3-9+tm1+(1:tm2),i1)=flag;
                end
            end
    end