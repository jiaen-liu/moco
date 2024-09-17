% 20230621: Jiaen Liu, normalize b1 data by the covariance matrix to improve sensitivity estimation, especially useful for 10.5T
function [b1,cp,cp_im]=ste_sense_ext(par,sorted_ste)
    if ~isfield(sorted_ste,'k_pimg')
        error('*** Parallel imaging data was not provided! ***');
    end
    para=sorted_ste.para;
    nch=para.n_channels;
    mid_pimg=par.use_ext_ref(1);
    % get data and parameters
    if ~isfield(sorted_ste.k_pimg,['mid' num2str(mid_pimg)])
        error(['*** Parallel imaging data for MID '...
               num2str(mid_pimg) ...
              ' was not provided! ***']);
    end
    k_pimg=getfield(sorted_ste.k_pimg,['mid' num2str(mid_pimg)]);
    para_pimg=getfield(sorted_ste.para_pimg,['mid' num2str(mid_pimg)]);
    %
    sitmp=size(k_pimg);
    im_pimg=fftmr(k_pimg,-1,[1:para_pimg.dimen])*prod(sitmp(1:para_pimg.dimen));

    if field_true(par,'ref_phase_b0')
        if para_pimg.b_epi_positive
            b0_data=squeeze(sum(im_pimg.*conj(im_pimg(:,:,:,:,1)),4));
            te=para_pimg.te_contr(:)*1e-3;
        elseif length(para_pimg.te_contr)>2
            b0_data=squeeze(sum(im_pimg(:,:,:,:,1:2:end).*conj(im_pimg(:,:,:,:,1)),4));
            te=para_pimg.te_contr(1:2:end)*1e-3;
        elseif length(para_pimg.te_contr)==2
            % correction of phase difference
            % only in the readout direction
            if abs(para_pimg.frequency/42.58e6-7)<0.5 && ...
                    isfield(par,'vendor') && ...
                    strcmp(par.vendor,'philips')
                idx=[1:8:nch]+[0:3].';
                idx1=idx(:);
                idx=[5:8:nch]+[0:3].';
                idx2=idx(:);
                idx=[idx1,idx2];
            else
                idx=[1:nch].';
            end
            prof=mean(mean(im_pimg,2),3);
            prof=squeeze(prof);
            prof=permute(prof,[1,3,2]);
            for igrp=1:size(idx,2)
                pd=dpOddEven(prof(:,:,idx(:,igrp)),2);
                im_pimg(:,:,:,idx(:,igrp),1:2:end)=...
                    im_pimg(:,:,:,idx(:,igrp),1:2:end).*...
                    exp(-1i*pd/2);
                im_pimg(:,:,:,idx(:,igrp),2:2:end)=...
                    im_pimg(:,:,:,idx(:,igrp),2:2:end).*...
                    exp(1i*pd/2);
            end
            b0_data=squeeze(sum(im_pimg.*conj(im_pimg(:,:,:,:,1)),4));
            te=para_pimg.te_contr*1e-3;
        end
        b0=b0_map(b0_data,te);
        ref_pha=exp(1i*b0*2*pi*te(1));
    else
        ref_pha=[];
    end
    
    im_pimg=squeeze(im_pimg(:,:,:,:,1,1,1,1));
    
    imc_pimg=sum(abs(im_pimg),4);
    % 0.3 and 0.01 for 10.5 T
    if isfield(par,'ref_mask_ratio')
        ref_mask_ratio=par.ref_mask_ratio;
    else
        ref_mask_ratio=0.1;
    end
    [~,thrd]=mask1d(imc_pimg(:),0.3,ref_mask_ratio);
    % [~,thrd]=mask1d(imc_pimg(:),0.3,0.01);
    % mask_pimg=imc_pimg>par.use_ext_ref(2);
    mask_pimg=imc_pimg>thrd;
    if field_true(sorted_ste,'sense_philips')
        im_pimg=permute(im_pimg,[2,3,1,4]);
        mask_pimg=permute(mask_pimg,[2,3,1]);
    end
    if ~isfield(par,'frac_sense_sm')
        par.frac_sense_sm=8;
    end
    if ~isfield(par,'cp_z_shift')
        par.cp_z_shift=0;
    end
    [b1,ref,cp,~,cp_im]=sense_m(im_pimg,mask_pimg,eye(nch,nch),...
                                'ress',para_pimg.ress,...
                                'frac_sense_sm',par.frac_sense_sm,...
                                'cp_z_shift',par.cp_z_shift,...
                                'ref_pha_in',ref_pha,...
                                'poly_fit',par.poly_fit,...
                                'poly_fit_ker',par.poly_fit_ker);
    
    if field_true(sorted_ste,'sense_philips')
       b1=permute(b1,[3,1,2,4]);
    end
    b1=covNorm(b1,inv(sorted_ste.inv_cov),4);
% $$$     if para.freq/42.58e6>6
% $$$         
% $$$     end
end