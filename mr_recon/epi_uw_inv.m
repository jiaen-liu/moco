function im_recon=epi_uw_inv(kd,para,b0,non_iter)
% mask_te_kps nrxnpxnsxn_te
    i2pi=1i*2*pi;
    
    nr=para.nr;
    np=para.np;
    ns=para.n_partitions;
    rph=para.sense_rate_p;
    rs=para.sense_rate_s;
    sshift=para.sshift;
    dph=para.dph;
    ds=para.ds;
    n=nr*np*ns;
    nps=np*ns;
    ninterl=para.n_interleaves;
    np_shot=np/ninterl/rph;
    nte=np_shot;
    nch=para.n_channels;
    es=para.echo_spacing*1e-6; % second
    te=para.te*1e-3+([0:nte-1]-floor(nte/2))*es;
    npk=np/rph;
    nsk=ns/rs;
    nkps=npk*nsk;
    nkps_ch=nkps*nch;
    %
    kd=reshape(kd,[nr,npk,nsk,nch]);
    % 
    sensit=para.b1;
    sensit=reshape(sensit,[nr,np,ns,nch]);
% $$$     sensit=reshape(sensit,[nr*np*ns,nch]);
% $$$     sensit=sensit*conj(chol(inv_cov,'lower'));
% $$$     sensit=sensit/prctile(abs(sensit(:)),95);

    %%
    % regularize sensit maps
    % abs_sensit=abs(sensit);
    % max_sensit=max(abs_sensit(:));    
    % sens_mask=abs_sensit<(max_sensit/100);
    % sensit(sens_mask)=max_sensit/100;
    % sens_mask=sum(sens_mask,4)==nch;
    %%
    % phasor 
    p=exp(i2pi*...
          (reshape(b0,[nr,np,ns]).*...
           reshape(te,[1,1,1,nte])));

    % For each TE, calculate the matrix of
    % FMk=Fh*Mk^h*Mk*F nxn
    % FM=sum over k(pk^h*FMk*pk)
    % pk is exp(i*2*pi*te_k*B0)
    % the calculation is compatible with caipirina
    % need a few parameters
    % rph,rs, sshift, dph, ds, n_interleaves
    % FM is stored as sparse matrix
    
    % np_shot
    y=[0:np-1]-floor(np/2);
    [Y1,Y2]=ndgrid(y,y);
    YmY=Y1-Y2;
    FMk=cell(nte,1);
    kp_te=zeros(ninterl,nte);
    ks_te=zeros(ns/rs,nte);
    %
    kd_zero=zeros(nr,np,ns,nch,nte);
    %
    for k=1:nte
        % get the kp positions for this echo
        kp_te(:,k)=[0:ninterl-1]*rph+(k-1)*ninterl*rph-floor(np/2)+dph;
        sshift_te=ds+mod((k-1)*sshift,rs);
        ks_te(:,k)=[0:ns/rs-1]*rs-floor(ns/2)+sshift_te;
        FMkp=sum(exp((-i2pi/np)*(YmY.*reshape(kp_te(:,k),[1,1,ninterl]))),3);
        FMk{k}=sparse(nps,nps);
        for i=1:ns
            for j=1:rs
                itmp=mod(i+ns/rs*(j-1)-1,ns);
                FMk{k}((i-1)*np+1:i*np,...
                       itmp*np+1:(itmp+1)*np)=FMkp*...
                    exp(-i2pi*floor(ns/2)*(j-1)/rs)*...
                    exp(-i2pi*sshift_te*(j-1)/rs);
            end
        end
        % zero fill the k space not acquired at the current TE
        kd_zero(:,mod(kp_te(:,k)+floor(np/2),np)+1,...
                mod(ks_te(:,k)+floor(ns/2),ns)+1,:,k)=...
            kd(:,mod(kp_te(:,k)+floor(np/2),np)+1,...
               mod(ks_te(:,k)+floor(ns/2),ns)+1,:);
    end
    imd_zero=fftmr(kd_zero,-1,[1,2,3])*(nr*np*ns);
    clear kd_zero;
    b=sum(conj(sensit).*...
          sum(reshape(conj(p),[nr,np,ns,1,nte]).*imd_zero,5),4);
    b=double(b);
    % need to choose full fov or accelerated imaging.
    AA=sparse(n,n);
    for ir=1:nr
        pr=squeeze(p(ir,:,:,:));
        PFM=sparse(nps,nps);
        pr_sp=cell(nte,1);
        for k=1:nte
            pr_sp{k}=spdiags(a2v(p(ir,:,:,k)),0,nps,nps);
        end
        sr_sp=cell(nch,1);
        for ich=1:nch
            sr_sp{ich}=spdiags(a2v(sensit(ir,:,:,ich)),0,nps,nps);
        end
        for k=1:nte
            PFM=PFM+conj(pr_sp{k})*FMk{k}*pr_sp{k};
        end
        PFM=PFM/double(ninterl*nte);
        A=sparse(nps,nps);
        for ich=1:nch
            A=A+conj(sr_sp{ich})*PFM*sr_sp{ich};
        end
        AA((ir-1)*nps+1:ir*nps,(ir-1)*nps+1:ir*nps)=A;
    end
    b=permute(b,[2,3,1]);
    im_recon=zeros(np*ns*nr,1);
    if ~non_iter
        [im_recon,flag]=pcg(AA,a2v(b),1e-3,3,[],[],zeros(n,1));
        [im_recon,flag]=pcg(AA,a2v(b),1e-3,3,[],[],im_recon);
    else
        % directly solving inverse problem
        % using a thresholded regularization
        v=spdiags(AA,0);
        vn=v*1.5;
        rr=50;
        thrd=max(abs(v))/rr;
        m=find(abs(v)<thrd);
        vn(m)=thrd;
        if thrd>0
            im_recon=(AA+spdiags(vn-v,0,n,n))\a2v(b);
        end
        % im_recon=a2v(b)./vn;
% $$$         tmp=v(find(v~=0));
% $$$         lambda=0.8;
% $$$         im_recon=(AA+lambda*speye(n,n))\a2v(b);
    end
    im_recon=reshape(im_recon,[np,ns,nr]);
    im_recon=permute(im_recon,[3,1,2]);
end
