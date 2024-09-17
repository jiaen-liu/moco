function [y,para,l,cov_mat]=recon_bmir_epi(mid,varargin)
    p=inputParser;
    apodiz=0.25;
    k_return=0;
% $$$     no_comb=0;
    no_pc=0;
    no_fov_crct=0;
    ord_pha_crct=2;
% $$$     mid_b1=[];
    mix=0;
    addParameter(p,'apodiz',apodiz,@isnumeric);
    addParameter(p,'k_return',k_return,@isnumeric);
% $$$     addParameter(p,'no_comb',no_comb,@isnumeric);
    addParameter(p,'no_pc',no_pc,@isnumeric);
    addParameter(p,'no_fov_crct',no_fov_crct,@isnumeric);
    addParameter(p,'ord_pha_crct',ord_pha_crct,@isnumeric);
    addParameter(p,'mix',mix,@isnumeric);
% $$$     addParameter(p,'mid_b1',mid_b1,@isnumeric);
    p.parse(varargin{:});
    apodiz=p.Results.apodiz;
    k_return=p.Results.k_return;
% $$$     no_comb=p.Results.no_comb;
    no_pc=p.Results.no_pc;
    no_fov_crct=p.Results.no_fov_crct;
    ord_pha_crct=p.Results.ord_pha_crct;
    mix=p.Results.mix;
% $$$     mid_b1=p.Results.mid_b1;
    if k_return
        no_comb=1;
    end
    %% set parameters
    para=extract_para_philips(mid);
    nr=para.nr;
    nr_os=para.nr_os;
    np=para.np;
    nch=para.n_channels;
    necho=para.nte_contr;
    nk_shot=para.nk_shot;
    n_tot_line_shot=nk_shot*necho;
    ns=para.n_slices;
    npar=para.n_partitions;
    s1=para.sense_rate_p;
    s2=para.sense_rate_s;
    ninterl=para.n_interleaves;
    nps=np/s1;
    npars=npar/s2;
    nreps=para.n_reps;
    
    tfe_factor=para.tfe_factor;
    n_delays=para.n_delays;
    % get noise data
    n=read_raw_philips(mid,'type',5);
    cov_mat=cov(conj(n));
    cov_mat=cov_mat/2;
    % get image data
    [dx,l]=read_raw_philips(mid,'type',1,'corr_half_fov',1,...
                            'mix',mix);
    if para.k_ellip
        npe=length(l.ky(:))/necho;
    else
        npe=nps*npars;
    end
    
    % regridding
    if para.ramp_samp_frac > 0.01
        dx=calc_regrid_mat(para.nus_enc,nr_os,apodiz)*dx;
    else
        dx=apodize_arr(dx,apodiz,1);
        dx=fftmr(dx,-1,1);
    end
    % truncate readout oversampling
    dx=dx(idx_truncate(size(dx,1),nr),:);
    % correct fov in phase encoding direction
    % note: fov is not corrected for 3d navigator in all three directions
    if ~no_fov_crct
        if ~para.isgre
            % gre data has been corrected by the scanner
            % no correction is needed
            % epi is only corrected in the slice direction
            dx=reshape(dx,[nr,nch,numel(dx)/nr/nch]);
            dx=dx.*...
               exp(-1i*reshape(2*pi*l.ky*para.p_shift(1)/para.sampled_fovs(2),...
                               [1,1,numel(dx)/nr/nch]));
        end
    end
    if (ns>1 && tfe_factor>1) || (ns>1 && n_delays>1)
        error('*** Multi-slice TFE or multi-delay not supported! ***');
    end
    % reshape dx
    n_tr=numel(dx)/nr/nch/n_tot_line_shot/ns;
    dx=reshape(dx,[nr,nch,n_tot_line_shot,ns,n_tr]);
    % apodize k-space
    if nps>2
        dx=dx.*reshape(win_tukey(nps,apodiz,(l.ky-min(l.ky(:)))/nps),...
                       [1,1,n_tot_line_shot,ns,n_tr]);
    end
    if npars>2
        dx=dx.*reshape(win_tukey(npars,apodiz,(l.kz-min(l.kz(:)))/npars),...
                       [1,1,n_tot_line_shot,ns,n_tr]);
    end
    % correct slice order
    loca=reshape(l.loca,[n_tot_line_shot,ns,n_tr]);
    loca=reshape(loca(1,:,:),[ns,n_tr]);
    dxcp=dx;
    for i=1:n_tr
        dx(:,:,:,loca(:,i)+1,i)=dxcp(:,:,:,:,i);
    end
    clear dxcp;
    % nr x n_tot_line_shot x ns x nch x nshot
    dx=permute(dx,[1,3,4,2,5]);
    if ~para.isgre && abs(para.frequency/42.58e6-7)<0.5
        dx=reshape(dx,[nr,nk_shot,necho,ns,nch,n_tr]);
        dx(:,2:2:end,:,:,:,:)=...
            -dx(:,2:2:end,:,:,:,:);
        dx=reshape(dx,[nr,n_tot_line_shot,ns,nch,n_tr]);
    end
    % correct epi ghost
    if ~para.isgre && ~no_pc
        disp('*** EPI ghost correction ***');
        % retrieve blipoff data for epi ghost correction
        if abs(para.frequency/42.58e6-7)<0.5
            idx=[1:8:nch]+[0:3].';
            idx1=idx(:);
            idx=[5:8:nch]+[0:3].';
            idx2=idx(:);
            idx=[idx1,idx2];
        else
            idx=[1:nch].';
        end
        dblpo=get_phc_philips(mid);

        for i=1:size(idx,2)
            for ie=1:necho
                dx_contr=dx(:,(ie-1)*nk_shot+1:ie*nk_shot,:,idx(:,i),:);
                dblpo_contr=dblpo(:,(ie-1)*nk_shot+1:...
                                  ie*nk_shot,idx(:,i),:);
                dx_contr=pha_crct_epi(dx_contr,dblpo_contr,...
                                      para.ro_pol(:,ie),para.ro_pol(:,ie),...
                                      ord_pha_crct,1);
                dx(:,(ie-1)*nk_shot+1:ie*nk_shot,:,idx(:,i),:)=dx_contr;
            end
        end
    end
    % correct fov in readout direction
    % note: fov is not corrected for 3d navigator in all three directions
    if ~no_fov_crct
        if ~para.isgre
            % epi is only corrected in the slice direction
            lin_pha_ro=([0:nr-1].'-floor(nr/2))*2*pi*para.m_shift(1)/para.fovr;
            dx=fftmr(dx,1,1).*exp(-1i*lin_pha_ro);
            dx=fftmr(dx,-1,1);
        else
            % gre data has been corrected by the scanner
            % with a residual
            m_freq_factor_resi=-(para.m_freq_factor(1)-para.m_shift(1)*42.58);
            lin_pha_ro=([0:nr-1].'-floor(nr/2))*2*pi*m_freq_factor_resi/42.58/para.fovr;
            dx=fftmr(dx,1,1).*exp(-1i*lin_pha_ro);
            dx=fftmr(dx,-1,1);
        end
    end
    % return k-space data for later processing
    if k_return
        y=fftmr(dx,1,1);
        l.ky=l.ky*s1;
        l.kz=l.kz*s2;
        return;
    end
    % generate an uncombined image without noise normalization
    ky=reshape(l.ky,[nk_shot,necho,ns,n_tr]);
    kz=reshape(l.kz,[nk_shot,necho,ns,n_tr]);
    kymin=min(ky(:));
    kzmin=min(kz(:));
    y=zeros(nr,nch,ns,necho,nps*npars,nreps);
    dx=reshape(dx,[nr,nk_shot,necho,ns,nch,n_tr]);
    % nr x nch x ns x echo x npe x nreps
    dx=permute(dx,[1,5,4,3,2,6]);
    dx=reshape(dx,[nr,nch,ns,necho,npe,nreps]);
    ky=permute(ky,[3,2,1,4]);
    ky=reshape(ky,[ns,necho,npe,nreps]);
    kz=permute(kz,[3,2,1,4]);
    kz=reshape(kz,[ns,necho,npe,nreps]);
    % reorder phase encoded lines
    for irep=1:nreps
        for iecho=1:necho
            for is=1:ns
                ipe1=ky(is,iecho,:,irep)-kymin+1;
                ipe2=kz(is,iecho,:,irep)-kzmin+1;
                ipe=ipe1+(ipe2-1)*nps;
                ipe=ipe(:);
                y(:,:,is,iecho,ipe,irep)=...
                    dx(:,:,is,iecho,:,irep);
            end
        end
    end
    y=reshape(y,[nr,nch,ns,necho,nps,npars,nreps]);
    % nr, npe1, npe2, ns, nch, necho, nreps
    y=permute(y,[1,5,6,3,2,4,7]);
    y=fftmr(y,-1,[2,3]);
end