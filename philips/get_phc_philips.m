function ph=get_phc_philips(mid)
    apodiz=0.25;
    para=extract_para_philips(mid);
    nr_os=para.nr_os;
    nr=para.nr;
    nch=para.n_channels;
    necho=para.nte_contr;
    nk_shot=para.nk_shot;
    n_tot_line_shot=nk_shot*necho;
    ns=para.n_slices;
    [ph,l]=read_raw_philips(mid,'type',3);
    nrep_ph=numel(unique(l.aver));
    tfe_factor=para.tfe_factor;
    n_delays=numel(unique(l.card));
    if para.ramp_samp_frac > 0.01
        ph=calc_regrid_mat(para.nus_enc,nr_os,apodiz)*ph;
    else
        ph=apodize_arr(ph,apodiz,1);
        ph=fftmr(ph,-1,1);
    end
    ph=ph(idx_truncate(size(ph,1),nr),:);
    if (ns>1 && tfe_factor>1) || (ns>1 && n_delays>1)
        error('Multi-slice TFE or multi-delay is not supported!');
    end
    ph=reshape(ph,[nr,nch,n_tot_line_shot,ns,tfe_factor,n_delays,nrep_ph]);
    % correct slice order
    loca=reshape(l.loca,[n_tot_line_shot,ns,tfe_factor,n_delays,nrep_ph]);
    loca=loca(1,:,1,1,1);
    ph_cp=ph;
    ph(:,:,:,loca+1,:,:,:)=ph_cp;
    clear ph_cp;
    % nr x nline x nch x tfe_factor x ndelays x ns x rep
    ph=permute(ph,[1,3,2,5,6,4,7]);
    ph=reshape(ph,[nr,n_tot_line_shot,nch,tfe_factor*n_delays*ns*nrep_ph]);
    if ~para.isgre && abs(para.frequency/42.58e6-7)<0.5
        ph(:,2:2:end,:,:)=...
            -ph(:,2:2:end,:,:);
    end
end