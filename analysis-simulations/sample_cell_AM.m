clear all
tic
n_simu=1;
i0=1;
nc=2e4;% number of cells mimiced
nic=21;% number of initial values for chromatin accessibility (nucleosome opening) (0:0.05:1)*9
nI=20; % number of light intensities scanned
np=15; % number of values for HDC ?HDAC4/5 molecules) or Gt (total GAVPO molecules) 
nv=9;  % number of variables in the model

out_total=zeros(nic, nI,  np, np, nv,n_simu);   
param_total=zeros(np*np,26,n_simu);
ic_total=zeros(nI,nic,nv,np*np,n_simu);
fileinfo1=fileinfo(j00);
    filename='AMV40-simu-result.mat';
    load(filename);
    I=result.I;nI=length(I);
    ic=result.ic;
    out=result.outsIC;
    out2=permute(out(:,:,end,:,:),[1 2 5 4 3 ]);
    out21=reshape(out2,nic,nI,np,np,nv);
    param=result.param;

    out_total(:,:,:,:,:,i0)=out21;
    param_total(:,:,i0)=param;
    ic_total(:,:,:,:,i0)=ic;
    filemat='AM-simu-mimic-cell-populations.mat';
    [I]=minic_cell_population_sampling(filemat,i0,out_total,param_total,ic_total,nc,nic,nI,np,nv,I);
