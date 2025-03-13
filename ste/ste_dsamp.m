function [para,kd,par]=ste_dsamp(para,kd,par)
    mask_empty_cl=false(par.nc,1);
    ncontr=length(par.icontr);
    npe=numel(kd)/para.nr/para.n_channels/ncontr;
    sik=size(kd);
    kd=reshape(kd,[para.nr,npe,para.n_channels,ncontr]);
    mask_kd=false(npe,1);
    nk=0;
    for i=1:par.nc
        mask_gre_cl=(mod(para.md(i).k(1,:)+floor(para.np/2),par.dsamp(1))==0) & ...
            (mod(para.md(i).k(2,:)+floor(para.n_partitions/2),par.dsamp(2))==0);

        mask_kd(nk+1:nk+para.md(i).nk)=mask_gre_cl;
        nk=nk+para.md(i).nk;
        if ~any(mask_gre_cl)
            mask_empty_cl(i)=true;
        else
            para.md(i).nk=total(mask_gre_cl);
            para.md(i).k=para.md(i).k(:,mask_gre_cl);
            para.md(i).m=para.md(i).m(:,:,mask_gre_cl);
            para.md(i).gb0=para.md(i).gb0(:,mask_gre_cl);
            para.md(i).db0=para.md(i).db0(mask_gre_cl);
        end
    end
    %% delete empty md
    para.md(mask_empty_cl)=[];
    para.nm=para.nm-total(mask_empty_cl);
    par.nc=par.nc-total(mask_empty_cl);
    %% downsample kd
    kd=kd(:,mask_kd,:,:);
    npe_ds=total(mask_kd);
    kd=reshape(kd,[para.nr*npe_ds,para.n_channels,ncontr]);
end