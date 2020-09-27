function [I]=minic_cell_population_sampling(filemat,i0,out_total,param_total,ic_total,nc,nic,nI,np,nv,I)
rng('shuffle');
x5=-5:0.1:7; nx5=length(x5);              

%simulated result, 
% 21 initial nucleosome conditions,nic
% 20 light intensities, nI
% 15 initial values for HDAC, np
% 15 initial values for total GAVPO. np
% 9 vairiables, nv
out21=out_total(:,:,:,:,:,i0);
param=param_total(:,:,i0);
ic=ic_total(:,:,:,:,i0);

CO_b=0.1; % poisson sampling of initial nucleosome openning for individual cells, sigma=CO_b*(nic-1)=6;
m_H=8;  % the median from the simulations with 15 HDAC values
s_H=0.05; % sigma for sampling HDC values for individual cells
m_G=6;  % the median from the simulations with 15 GAVPO values
s_H=0.15; % sigma for sampling total GAVPO values for individual cells

    % initial values for nucleosome
    % 9 sites on genome, each has 2, so max=18
    [ics0]=poisson_nucleosome(nc,nic,CO_b);
    
    % initial values for HDAC
    [ics2]=normal_HDAC_offset(nc,s_H,15,m_H,7);
    
    % initial value for GAVPO
    [ics3]=normal_GAVPO(nc,s_G,15,m_G,14);
    
    ic=(ics3-1)*nic*np+(ics2-1)*nic+ics0;
 
    % out24n, mimic cell populations (nc=2e4) from simulated data
    out22=permute(out21,[1 3 4 2 5]);
    out23=reshape(out22,nic*np*np,nI,nv);
    out24n=out23(ic,:,:);                    
    for i1=1:9
        % add log noise, 1% of the max values for each variables, cv=0.4
        mx=mean(out24n(:,end,i1));
        addn=lognrnd(1,1.2,[nc,20]); % cv=0.4
        addn=addn/mean(addn);
        out24n(:,:,i1)=out24n(:,:,i1)+addn*mx*0.01;
        
        % multiply log noise, cv=0.4
        addn=lognrnd(1,1.2,[nc,20]);% cv=0.4
        addn=addn/mean(addn);
        out24n(:,:,i1)=out24n(:,:,i1).*addn;
    end
    % histograms for all 9 variables for all 11 light intensities
    hgn=zeros(nx5,nI,nv);
    for i1=1:nI
        for i2=nv
            hgn(:,i1,i2)=hist(log10(out24n(:,i1,i2)),x5)/nc;
        end
    end

save(filemat,'out24n','hgn','param','ic','x5');
end

function [ics0]=poisson_nucleosome(nc,nic,delta)
    ics=poissrnd(delta*(nic-1),[10*nc,1])/((nic-1)/2);
    ics1=ceil(ics*nic/2);
    j1=find(ics1>nic);
    ics1(j1)=[];
    j2=find(ics1<1);
    ics1(j2)=[];
    ics0=(ics1(1:nc));
end
function [ics2]=normal_HDAC_offset(nc,sig,np,npc,npw)
    ics=randn(nc*10,1)*sig+1;
    ics1=ceil(ics*npw/2+(npc-npw/2-0.5));
    j1=find(ics1>np);
    ics1(j1)=[];
    j2=find(ics1<1);
    ics1(j2)=[];
    ics2=(ics1(1:nc));
end
function [ics3]=normal_GAVPO(nc,sig,np,npc,npw)
    ics=randn(nc*10,1)*sig+1;
    ics1=ceil(ics*npw/2+(npc-npw/2-0.5));
    j1=find(ics1>np);
    ics1(j1)=[];
    j2=find(ics1<1);
    ics1(j2)=[];
    ics3=(ics1(1:nc));
end