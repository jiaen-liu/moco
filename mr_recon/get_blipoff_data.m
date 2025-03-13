function dblpo=get_blipoff_data(mid)
apodiz=0.25;
para=extract_para(mid);
ro_pol=para.ro_pol_ref(:);
nr_os=para.nr_os;
nr=para.nr;
n_echo_contr=para.np/para.n_interleaves/...
    para.sense_rate_p;
n_echo_shot=0;
ncontr=para.n_contrasts;
nch=para.n_channels;
npar=para.n_partitions;
s2=para.sense_rate_s;
ns=para.n_slices;
for i=1:ncontr
    n_echo_shot=n_echo_shot+n_echo_contr+...
        para.n_refs(i);
end
dblpo=readSortSiem(mid,'blipoff',1);
% correct polarity
dblpo(:,find(ro_pol==0),:,:,:)=...
    flipdim(dblpo(:,...
                  find(ro_pol==0),...
                  :,:,:),1);
npar_blpo=size(dblpo,5);
% regridding
dblpo=regridding_arr(dblpo,para,apodiz);
dblpo=dblpo(idx_truncate(nr_os,nr),:,:,:,:);

dblpo=reshape(dblpo,[nr,n_echo_shot,ns,...
                    nch,npar_blpo]);
% conver to image domain
dblpo=fftmr(dblpo,-1,5);
dblpo=permute(dblpo,[1,2,4,5,3]);
% Based on the obersevation that 
% eddy current is mostly along read-out direction
% other dimentions will be averaged in dpOddEven.m
dblpo=mean(dblpo,6);
dblpo=reshape(dblpo,[nr,n_echo_shot,...
                    nch,ns*npar_blpo]);