function [im_uw,b0,b1_o,b1_o_nm]=ste_uw_inv_mv(im_input,n_interl,te,...
                                                  mask_thrd,...
                                                  pe_dir,...
                                                  cov_mat,...
                                                  niter,...
                                                  gre_calib,...
                                                  gre_b0,...
                                                  b1_input,...
                                                  non_iter)

    % May need to iterate a few steps
    nch=size(cov_mat,1);
    kd=fftmr(im_input,1,[1,2,3]);
    if nch>1
        [nr,np,ns,nch,ncontr,nv]=size(im_input);
    else
        [nr,np,ns,ncontr,nv]=size(im_input);
    end
    im_input=reshape(im_input,[nr,np,ns,nch,ncontr,nv]);
    n=nr*np*ns;
    nline_interl=np/n_interl;
    im_uw=zeros(nr,np,ns,ncontr,nv);
    inv_cov=inv(cov_mat);
    % extract sensitivity
    % save the data and let idl process it
    if numel(mask_thrd)>1
        mask=mask_thrd;
    else
        if isempty(gre_calib)
            mask=squeeze(sum(abs(im_input(:,:,:,:,1,:)),4))>mask_thrd;
        else
            mask=sum(abs(gre_calib),4)>mask_thrd;
        end
    end

    
    if isempty(b1_input)
        
% $$$         fn_sen=rp('sens_ste4idl.svd','test');
% $$$         fn_out=rp('sens_ste4matlab.svd','test');
% $$$         idl_e='$IDL_DIR/bin/idl -quiet -e ';
% $$$         if isempty(gre_calib)
% $$$             tmp=struct('mask',mask,'im',squeeze(im_input(:,:,:,:,1,:)),'fn_out',fn_out);
% $$$         else
% $$$             tmp=struct('mask',mask,'im',gre_calib,'fn_out',fn_out);
% $$$         end
% $$$         save_data(fn_sen,tmp);
% $$$         cmd=[idl_e '''' 'sense4matlab,' '"' fn_sen '"' ''''];
% $$$         disp(cmd);
% $$$         status=system(cmd);
% $$$         if status ~= 0
% $$$             error('Error in calculating sensitivity maps using IDL!');
% $$$         end
% $$$         s=read_data(fn_out);
% $$$         b1=s.sensit;        
        [b1,ref]=sense_m(squeeze(im_input(:,:,:,:,1,:)),...
                         mask,...
                         eye(nch,nch),'cp_all',1);

        b1_o_nm=b1;
        % prepare input to recon_mb_gre
        b1=reshape(b1,[nr*np*ns,nch,nv]);
        b1=permute(b1,[1,3,2]);
        b1=reshape(b1,[nr*np*ns*nv,nch]);
        b1=b1*conj(chol(inv_cov,'lower'));
        b1=reshape(b1,[nr,np,ns,nv,nch]);
        b1=permute(b1,[1,2,3,5,4]);
        b1=b1/prctile(abs(b1(:)),95);
    else
        b1=b1_input;
        b1_o_nm=[];
    end
    b1=reshape(b1,[nr,np,ns,nch,nv]);
    b1_o=b1;
    if isempty(niter)
        niter=0;
    end
    % extract initial b0
    im_cmb=squeeze(sum(im_input.*conj(im_input(:,:,:,:,1,:)),4));
    if pe_dir(1)~=pe_dir(2)
        dte=(te(end/2)-te(1))*1e-3;
    else
        dte=(te(end/2+1)-te(1))*1e-3;
    end
    if isempty(gre_b0)
        ph=squeeze(angle(im_cmb(:,:,:,end,:)));
        cp=squeeze(mean(reshape(ph(floor(nr/2):floor(nr/2)+2,...
                                   floor(np/2):floor(np/2)+2,...
                                   floor(ns/2)+1,:),[9,nv]),1));
        ph=ph-reshape((cp-cp(1)-wrapToPi(cp-cp(1))),[1,1,1,nv]);

        b0=zeros(nr,np,ns,nv);
        for i=1:nv
            b0(:,:,:,i)=unwrapper_3d_mask(ph(:,:,:,i))/...
                2/pi/dte;
        end
        b0=b0.*mask;
    else
        cp=0;
        b0=gre_b0;
    end
    sik=size(kd);
    kd=reshape(kd,[n,nch,ncontr,nv]);
    kd=permute(kd,[1,3,4,2]);
    kd=reshape(kd,[n*ncontr*nv,nch]);
    kd=kd*conj(chol(inv_cov,'lower'));
    kd=reshape(kd,[n,ncontr,nv,nch]);
    kd=permute(kd,[1,4,2,3]);
    im_uw=zeros(nr,np,ns,ncontr,nv);
    nv_seg=40;
    nseg=ceil(nv/nv_seg);
    idx_seg=zeros(2,nseg);
    for iseg=1:nseg
        idx_seg(1,iseg)=nv_seg*(iseg-1)+1;
        idx_seg(2,iseg)=min(nv_seg*iseg,nv);
    end
    
    kd_cell=cell(nseg);
    b1_cell=cell(nseg);
    for iseg=1:nseg
        kd_cell{iseg}=...
            kd(:,:,:,...
               idx_seg(1,iseg):idx_seg(2,iseg));
        b1_cell{iseg}=b1(:,:,:,:,idx_seg(1,iseg):idx_seg(2,iseg));
    end
    
    clear kd b1
    for iseg=1:nseg
        print_countdown(nseg,iseg,'*** B0 Unwarping: ');
        nv_seg_cur=idx_seg(2,iseg)-idx_seg(1,iseg)+1;
        kd_seg=kd_cell{iseg};
        b1_seg=b1_cell{iseg};
        b0_seg=b0(:,:,:,idx_seg(1,iseg):idx_seg(2,iseg));
        im_uw_seg=zeros(nr,np,ns,ncontr,nv_seg_cur);
        parfor iv=1:nv_seg_cur
            b0v=b0_seg(:,:,:,iv);
            para_recon=[];
            im_uw_v=zeros(nr,np,ns,ncontr);
            for iiter=1:niter
                for icontr=1:ncontr
                    te_contr=te((icontr-1)*nline_interl+1:icontr*nline_interl);
                    if pe_dir(icontr)==-1
                        te_contr=flip(te_contr);
                    end

                    kd_contr=kd_seg(:,:,icontr,iv);
                    kd_contr=reshape(kd_contr,[nr,np,ns,nch]);
                    
                    para_recon.nr=nr;
                    para_recon.np=np;
                    para_recon.n_partitions=ns;
                    para_recon.sense_rate_p=1;
                    para_recon.sense_rate_s=1;
                    para_recon.sshift=0;
                    para_recon.dph=0;
                    para_recon.ds=0;
                    para_recon.n_channels=nch;
                    para_recon.echo_spacing=(te(2)-te(1))*1000;
                    para_recon.te=te_contr(end/2+1);
                    para_recon.n_interleaves=n_interl;
                    para_recon.b1=b1_seg(:,:,:,:,iv);
                    im_uw_v(:,:,:,icontr)=epi_uw_inv(kd_contr,para_recon,b0v,non_iter);
                end
                ph_resi=angle(im_uw_v(:,:,:,2)./im_uw_v(:,:,:,1));
                b0_residual=unwrapper_3d_mask(ph_resi)/...
                    2/pi/dte;
                
                b0v=b0v+b0_residual;

                ph_cur=2*pi*b0v*dte;
                cp_cur=mean(a2v(ph_cur(floor(nr/2):floor(nr/2)+2,...
                                       floor(np/2):floor(np/2)+2,...
                                       floor(ns/2)+1)));
                ph_cur=ph_cur-(cp_cur-cp(1)-wrapToPi(cp_cur-cp(1)));
                b0v=ph_cur/2/pi/dte;
                im_uw_seg(:,:,:,:,iv)=im_uw_v;
            end
            b0_seg(:,:,:,iv)=b0v;

        end
        b0(:,:,:,idx_seg(1,iseg):idx_seg(2,iseg))=b0_seg;
        im_uw(:,:,:,:,idx_seg(1,iseg):idx_seg(2,iseg))=im_uw_seg;
    end
end