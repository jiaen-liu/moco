% 2022-12-23: Jiaen Liu, added support for Philips
% 2023-06-21: Jiaen Liu, improve sensitiity map estimation for 10.5T
function [im_recon,par,para_mr]=recon_epi_ste_beta(mid,par_script,chain,varargin)

% test if INTEL MKL works
    if ~valid_mkl_sparse_mult
        error('*** Intel MKL was not verified! ***');
    else
        disp('*** Intel MKL was verified! ***');
    end
    % 
    par=ste_default_conf;
    if ischar(par_script)
        run(par_script);
    elseif isstruct(par_script)
        par=pass_var_struct(par,par_script);
    end
    if par.n_seg_mem>1 && par.nc_max>1
        error('*** When k-means is used, memory segmention is not allowed! ***');
    end
    % 
    if ~exist('chain')
        % chain: 'sense', 'motrec', 'frec', 'motest', 'fmotest', 'db0','cluster','grerec'
        % when number of cluster is changed, need to run 'motest',
        % 'frec', 'grerec'
        chain=[];
    elseif ismember('full',chain)
        chain={'sense', 'frec','motrec', 'motest', 'fmotest','db0','cluster','grerec'};
    end
    
    p=inputParser;
    sensit='reg';
    mot_uw=[];
    max_ord=5;
    addParameter(p,'sensit',sensit,@ischar);
    addParameter(p,'mot_uw',mot_uw,@isnumeric);
    addParameter(p,'max_ord',max_ord,@isnumeric);
    p.parse(varargin{:});
    sensit=p.Results.sensit;
    mot_uw=p.Results.mot_uw;
    max_ord=p.Results.max_ord;
    %
    im_recon=[];
    % ----------------------------------------------- %
    % create subfolder to store temporary files
    % ----------------------------------------------- %
    % subdir=rp('ste' int2str(mid)]);
    subdir=rp('');
    if exist(subdir) ~= 7
        mkdir(subdir);
    end
    % ----------------------------------------------- %
    % get the pre-processed ste data
    % ----------------------------------------------- %
    disp('*** Loading the pre-processed ste data ***');
    fn=rp(['mid' int2str(mid) '.steref4recon.svd']);
    data=read_data(fn);
    if ~isfield(data,'vendor')
        vendor='siemens';
    else
        vendor=data.vendor;
        if iscell(vendor)
            vendor=vendor{1};
        end
    end
    par.vendor=vendor;
    para_mr=data.para;
    if isfield(para_mr,'freq') && abs(para_mr.freq/42.53e6-10.5)<0.2
        par.frac_sense_sm=2;
    elseif isfield(para_mr,'freq') && para_mr.freq/42.53e6<7.5
        par.frac_sense_sm=8;
    else
        par.frac_sense_sm=8;
    end
    data=struct2double(data);
    % sort the ste data
    sorted_ste=sort_ste(data);
    nshot_fv=sorted_ste.nshot_fv;
    ncontr_ste=sorted_ste.ncontr;
    pe_dir=ste_pe_dir(data);
    nch=sorted_ste.para.n_channels;
    %% calculate certain parameters
    if sorted_ste.para.isgre
        if sorted_ste.para.n_contrasts>1
            necho=length(sorted_ste.para.te_contr);
        else
            necho=sorted_ste.para.n_refs(1)+1;
        end
        
    else
        necho=length(sorted_ste.para.te_contr);
    end
    npar=sorted_ste.para.n_partitions;
    npar_acc=npar/sorted_ste.para.sense_rate_s;
    n_interl=sorted_ste.para.n_interleaves;
    nreps=sorted_ste.para.n_reps;
    res_ste=[sorted_ste.para.steref_res_r,...
             sorted_ste.para.steref_res_p,...
             sorted_ste.para.steref_res_s];
    par.ste_intp_res=min(res_ste);
    nxste=sorted_ste.para.steref_dim_r;
    nyste=sorted_ste.para.steref_dim_p;
    nzste=sorted_ste.para.steref_dim_s;
    
    pix_ste_shift=[sorted_ste.para.steref_dim_r,...
                   sorted_ste.para.steref_dim_p,...
                   sorted_ste.para.steref_dim_s];
    % pix_ste_shift=pix_ste_shift-2*floor(pix_ste_shift/2)-1;
    pix_ste_shift=-mod(1+pix_ste_shift,2)/2;
    % temporary test
    % pix_ste_shift=[0,0,0];
    %
    if ~isfield(par,'icontr')
        par.icontr=[1:necho];
    else
        par.icontr=min(par.icontr,necho);
    end
    % use external or internal sensitivity reference
    if isfield(par,'use_ext_ref') 
        if par.use_ext_ref(1)==0
            par.use_ext_ref=[];
        elseif par.use_ext_ref(1)==1
            % use the 2d gre as reference
            par.use_ext_ref(1)=getMIDPI(data);
        end
    end
    mecho=par.mecho;
    % if number of ste contrast is 1, don't unwarp ste or
    % use multiple echo to calculate B0 change
    if data.pe.ncontr==1
        par.uw=0;
        par.fast_uw=0;
        mecho=0;
    end
    par.mecho=mecho;
    % on the other hand, if only single echo is used
    % unwarpping should be disabled.
    if mecho==0
        par.uw=0;
        par.fast_uw=0;
    end
    % ----------------------------------------------- %
    % Define file names
    % ----------------------------------------------- %
    if field_true(par,'n_any_correction') || ...
            field_true(par,'nb0')
        % No correction even the global B0 correction
        fn_gre=fullfile(subdir,['mid' int2str(mid) '.k_nnav_uncomb.svd']);
        % check if data exist
        if ~exist(fn_gre)
            error(['*** The image data doesn''t exist! ***' newline ...
                   fn_gre]);
        end
        if field_true(par,'n_any_correction')
            par.ncorrection=1;
        end
    else
        fn_gre=fullfile(subdir,['mid' int2str(mid) '.k_nav0_uncomb.svd']);
    end
    
    fn_final_prep=fullfile(subdir,['mid' int2str(mid) '.final_prep.svd']);
    
    if isfield(par,'use_ref_mid')
        fn_sen_1v=fullfile(subdir,['mid' int2str(par.use_ref_mid) '.sen_ste_1v.mat']);
    else
        fn_sen_1v=fullfile(subdir,['mid' int2str(mid) '.sen_ste_1v.mat']);
    end
% $$$     switch sensit
% $$$       case 'reg'
% $$$         fn_sen_1v=fullfile(subdir,['mid' int2str(mid) '.sen_ste_1v.svd']);
% $$$       case 'espi'
% $$$         fn_sen_1v=fullfile(subdir,['mid' int2str(mid) '.sen_ste_1v_espi.svd']);
% $$$     end
    fn_db0_cluster=fullfile(subdir,['mid' int2str(mid) '.db0_cluster.mat']);
    fn_b1=fullfile(subdir,['mid' int2str(mid) '.b1.mat']);
    if field_true(par,'cluster_b0')
        fn_db0_cluster=file_addext(fn_db0_cluster,'_cb0');
    else
        fn_db0_cluster=file_addext(fn_db0_cluster,'_cmt');
    end
    fn_db0_comb=fullfile(subdir,['mid' int2str(mid) '.db0_comb.mat']);
    fn_im_mot=fullfile(subdir,['mid' int2str(mid) '.mot_ste_rec.mat']);
    fn_im_mot_uw=file_addext(fn_im_mot,'_uw');
    if isfield(par,'use_ref_mid')
        fn_im_mot=file_addext(fn_im_mot,['_ref_mid' int2str(par.use_ref_mid)]);
    end
    
    if isfield(par,'use_ref_mid')
        fn_im_mot_uw=file_addext(fn_im_mot_uw,['_ref_mid' int2str(par.use_ref_mid)]);
    end
    fn_im_coreg_p=fullfile(subdir,['mid' int2str(mid) '.im_coreg_p.mat']);
    fn_im_coreg_f=fullfile(subdir,['mid' int2str(mid) '.im_coreg_f.mat']);
    if par.uw
        fn_im_coreg_p=file_addext(fn_im_coreg_p,'_uw');
        fn_im_coreg_f=file_addext(fn_im_coreg_f,'_uw');
    end
    fn_f=fullfile(subdir,['mid' int2str(mid) '.f_ste_rec.mat']);
    fn_f_uw=file_addext(fn_f,'_uw');
    if isfield(par,'use_ref_mid')
        fn_f=file_addext(fn_f,['_ref_mid' int2str(par.use_ref_mid)]);
    end
    
    if isfield(par,'use_ref_mid')
        fn_f_uw=file_addext(fn_f_uw,['_ref_mid' int2str(par.use_ref_mid)]);
    end
    fn_f_uncomb=fullfile(subdir,['mid' int2str(mid) '.f_ste_rec_uncomb.svd']);
    fn_f_uncomb_uw=file_addext(fn_f_uncomb,'_uw');
    fn_cp_ref=fullfile(subdir,['mid' int2str(mid) '.cp_ref.mat']);
% $$$     if par.mot_uw
% $$$         fn_f=file_addext(fn_f,'_uw');
% $$$     end
% $$$     if par.mot_uw
% $$$         fn_f_uncomb=file_addext(fn_f_uncomb,'_uw');
% $$$     end
    fn_motpar=fullfile(subdir,['mid' int2str(mid) '.mot_par.mat']);
    if field_true(par,'fast_uw')
        fn_motpar=file_addext(fn_motpar,'_uw');
    end
    fn_motpar_nuw=fullfile(subdir,['mid' int2str(mid) '.mot_par_nuw.mat']);
    fn_motpar_raw=fullfile(subdir,['mid' int2str(mid) '.mot_par_raw.mat']);
    fn_motpar_rf_uw=fullfile(subdir,['mid' int2str(mid) '.mot_par_rf_uw.mat']);
% $$$     if par.uw
% $$$         fn_motpar=file_addext(fn_motpar,'_uw');
% $$$     end
    fn_motparf=fullfile(subdir,['mid' int2str(mid) '.mot_par_f.mat']);
    fn_motparf_nuw=file_addext(fn_motparf,'_nuw');
    if par.uw
        fn_motparf=file_addext(fn_motparf,'_uw');
    end
    if isfield(par,'use_ref_mid')
        fn_grappa_calib=fullfile(subdir,['mid' int2str(par.use_ref_mid) '.ste_grappa_calib.mat']);
    else
        fn_grappa_calib=fullfile(subdir,['mid' int2str(mid) '.ste_grappa_calib.mat']);
    end
    
    fn_ste_info=fullfile(subdir,['mid' int2str(mid) '.ste_info.mat']);
    fn_motpar_comb=fullfile(subdir,['mid' int2str(mid) ...
                        '.mot_par_comb.mat']);
    if par.uw
        fn_motpar_comb=file_addext(fn_motpar_comb,'_uw');
    end
    fn_db0_f=fullfile(subdir,['mid' int2str(mid) '.db0_f.mat']);
    fn_db0_p=fullfile(subdir,['mid' int2str(mid) '.db0_p.mat']);
    if isfield(par,'use_ref_mid')
        fn_ste_mask=fullfile(subdir,['mid' int2str(par.use_ref_mid) '.ste_mask.mat']);
    else
        fn_ste_mask=fullfile(subdir,['mid' int2str(mid) '.ste_mask.mat']);
    end
    if mecho
        fn_db0_f=file_addext(fn_db0_f,'_me');
        fn_db0_p=file_addext(fn_db0_p,'_me');
    else
        fn_db0_f=file_addext(fn_db0_f,'_se');
        fn_db0_p=file_addext(fn_db0_p,'_se');
    end
    fn_b0_fit=fullfile(subdir,['mid' int2str(mid),'.db0_fit.mat']);
    fn_b0_fit_cluster=fullfile(subdir,['mid' int2str(mid),'.db0_fit_cluster.mat']);
    if field_true(par,'cluster_b0')
        fn_b0_fit=file_addext(fn_b0_fit,'_cb0');
    else
        fn_b0_fit=file_addext(fn_b0_fit,'_cmt');
    end
    fn_motpar_test=fullfile(subdir,['mid' int2str(mid),'.mot_par_test.mat']);
    if field_true(par,'cluster_b0')
        fn_motpar_test=file_addext(fn_motpar_test,'_cb0');
    else
        fn_motpar_test=file_addext(fn_motpar_test,'_cmt');
    end
    %%%%%%%%%%%%%%%%%5
    if field_true(par,'test_sense')
        % the data should be combination of static scans
        if ~sorted_ste.combined
            error('*** This is in the ''Test Sense'' mode. The data should be combined data! ***');
        end
        % only use full navigators
        par.mot_cmb=2; 

    end
    % check if 'isen' is defined
    if ~isfield(par,'isen')
        % fist check if it will be calculated in 'sense' step
        if ismember('sense',chain) && ~field_true(par,'auto_isen')
            error('*** A reference volume must be determined in ''sense''! ***');
        elseif ismember('sense',chain)
            disp('*** A reference volume will be determined in ''sense'' ***');
        elseif ~ismember('sense',chain)
            % load reference index from ste_info
            if ~exist(fn_ste_info)
                error('*** A reference volume must have been stored in a ste_info file! ***');
            end
            tmp=load(fn_ste_info);
            if isfield(tmp.ste_info,'isen')
                par.isen=tmp.ste_info.isen;
            else
                error('*** A reference volume must have been stored in ste_info if not to be determined in ''sense''! ***');
            end
        end

    end
    % check if full stes are acquired
    if field_true(par,'forceGRERef')
        fullSTEAcq=0;
    else
        fullSTEAcq=(sorted_ste.y_cyc~=0) & ...
            (sorted_ste.z_cyc~=0);
    end
    % for multi-repetition, only one cluster is allowed to
    % simplify workflow
    if sorted_ste.para.n_reps>1 && par.nc_max>1
        error('*** Only one cluster is allowed for repeated scans! ***');
    end
    % ----------------------------------------------- %
    % get a full-fov ste at certain time point
    % this volume will be the reference for the recon
    % ----------------------------------------------- %
    if ismember('sense',chain)
        % create markers for motion corrupted full aquisitions
        nvf_ceil=ceil(sorted_ste.nshot_nblpl/sorted_ste.nshot_fv);
        nvf=sorted_ste.nvf;
        fcrptd=zeros(nvf_ceil,1);
        freq_fluct=zeros(nvf_ceil,1);
        freq_fluct_strict=zeros(nvf_ceil,1);
        for i=1:nvf_ceil
            freq_blpl_fv_strict=sorted_ste.df_blpl(max((i-2)*sorted_ste.nblpl_fv+1,1):min((i+1)*sorted_ste.nblpl_fv+1,length(sorted_ste.df_blpl)));
            freq_blpl_fv=sorted_ste.df_blpl(max((i-1)*sorted_ste.nblpl_fv+1,1):min(i*sorted_ste.nblpl_fv+1,length(sorted_ste.df_blpl)));
            freq_fluct(i)=max(freq_blpl_fv)-min(freq_blpl_fv);
            freq_fluct_strict(i)=max(freq_blpl_fv_strict)-min(freq_blpl_fv_strict);
        end
        fcrptd=freq_fluct>par.max_f;
        if all(fcrptd)
            error('*** All full fovs are corrupted! ***');
        end
        if nvf_ceil>nvf
            fcrptd(nvf_ceil)=1;
        end
        ste_info.fcrptd=fcrptd;
        
        disp(['*** ' int2str(total(fcrptd)) '/' int2str(length(fcrptd)) ...
              ' volumes are corrupted by motion ***']);
        % determine reference volume
        if strcmp(vendor,'siemens')
            if field_true(par,'auto_isen')
                % find five full fovs with lowest freq_fluct
                par.isen=ste_isen(freq_fluct_strict,sorted_ste,par);
                disp(['*** Full ste volume ' int2str(par.isen) ' was choosen as the reference ***']);
                disp(['*** The frequency spans ' int2str(freq_fluct(par.isen)) ' Hz in this volume ***']);
            end
        elseif strcmp(vendor,'philips')
            % determine shots with central kyz
            [~,idx_kyz_0]=...
                min(col(min(sum(sorted_ste.kyz(:,1:sorted_ste.para.nk_shot,:).^2,1).^0.5,[],2)));
            par.isen=floor(idx_kyz_0/(sorted_ste.nshot_v+1))+1;
            par.isen=floor(par.isen/8)+1;
% $$$             par.isen=ste_isen(zeros(size(freq_fluct_strict)),sorted_ste,par);
        end
        ste_info.isen=par.isen;
        save_mat(fn_ste_info,'ste_info',...
                 'overwrite',1,'read4all',1);
        % Reference data for STE reconstruction
        if fullSTEAcq || strcmp(vendor,'philips')
            disp('*** Obtain a full fov for sensitivity/grappa-calib calculation ***');
            if strcmp(vendor,'philips')
                dk_sen=get_k_ste(sorted_ste,1,nshot_fv);
            else
                dk_sen=get_k_ste(sorted_ste,par.isen,nshot_fv);
            end
            dk_sen=dk_sen(:,:,:,:,1);
            sidk=size(dk_sen);
            im_sen=fftmr(dk_sen,-1,[1,2,3])*prod(sidk(1:3));
            if field_true(par,'shift_ste')
                im_sen=shift(im_sen,pix_ste_shift);
            end
        else
            % require external reference for
            % STE reconstruction
            [dk_sen,im_sen,para_sen]=ste_nav_ref(par,sorted_ste);
        end
        % generate mask
        if isfield(par,'use_saved_mask')
            if ~exist(par.use_saved_mask)
                error('*** External mask doesn''t exist! ***');
            end
            disp(['*** Loading mask from ',...
                  par.use_saved_mask,' ***']);
            load(par.use_saved_mask);
        else
            im_sen_mag=sum(abs(im_sen),4);
            % 0.3 and 0.03 for 10.5T data
            [~,mask_thrd]=mask1d(im_sen_mag(:),0.3,0.15);
            mask=im_sen_mag>mask_thrd;
            % extract brain mask
            mask_brain=extract_brain_mask(mask);            
        end
        [sensit,ref,cp]=sense_m(im_sen,mask,eye(nch,nch),'ress',sorted_ste.para.steref_res_s);
        save_mat(fn_sen_1v,'sensit','ref','cp',...
                 'overwrite',1,'read4all',1);
        nz=size(mask,3);
        mask_sh=mask;
        mask_sh(:,:,1:floor(end/4))=0;
        mask_sh(:,:,floor(end*3/4):end)=0;
        mask_brain_sh=mask_brain;
        mask_brain_sh(:,:,1:floor(nz/4))=0;
        mask_brain_sh(:,:,floor(nz*3/4):nz)=0;
% $$$         if ~isdeployed() && ~contains(getenv('HOSTNAME'),'biowulf','IgnoreCase',true)
% $$$             imshow3d(mask_brain);
% $$$         end
        save_mat(fn_ste_mask,'mask','mask_sh','mask_brain','mask_brain_sh',...
                 'overwrite',1,'read4all',1);
        % create calibration data for grappa reconstruction
        disp('*** Calibrating GRAPPA ***');
        if fullSTEAcq || strcmp(vendor,'philips')
            nxgc=nxste;
            nygc=nyste;
            nzgc=nzste;
        else
            nxgc=para_sen.nr;
            nygc=para_sen.np;
            nzgc=para_sen.n_partitions;
        end
        par_g.ncalx=nxgc;
        par_g.ncaly=nygc;
        par_g.ncalz=nzgc;
        par_g.nd=sorted_ste.para.dimen;
        par_g.ry=data.pe.sp;
        par_g.rz=data.pe.r/data.pe.sp;
        par_g.dz=data.pe.dkz;
        par_g.nx=nxste;
        par_g.ny=nyste;
        par_g.nz=nzste;
        if ~isfield(par,'nkx_g')
            par.nkx_g=3;
            par.nky_g=2;
            par.nkz_g=2;
        end
        par_g.nkx=par.nkx_g;
        par_g.nky=par.nky_g;
        par_g.nkz=par.nkz_g;
        [w,ker]=grappa_calib(dk_sen,par_g);
        save_mat(fn_grappa_calib,'w','ker','par_g',...
                 'overwrite',1,'read4all',1);
    end
    % ----------------------------------------------- %
    % Reconstruct full-fov images
    % ----------------------------------------------- %
    if ismember('frec',chain) && fullSTEAcq
        load(fn_ste_mask);
        load(fn_sen_1v);
        disp('*** Reconstructing full-fov images ***');
        imf=proc_sorted_ste(sorted_ste,'fullfov',1,'no_comb',1);
        imf(isnan(imf))=0;
        if field_true(par,'shift_ste')
            imf=shift(imf,pix_ste_shift);
        end
        imf_uncomb=imf;
        % save_data(fn_f_uncomb,imf_uncomb); % data is large, use Jacco's save_data
        % save_data(fn_f_uncomb_uw,'imf_uncomb_uw');
        [nx,ny,nz,nch,ncontr,nvf]=size(imf_uncomb);
        imf_mag=reshape(sum(abs(imf_uncomb),4),...
                        [nx,ny,nz,ncontr,nvf]);
        imf_pha=sum((abs(imf_uncomb(:,:,:,:,1,:)).*imf_uncomb).*...
                    exp(-1i*reshape(cp,[1,1,1,nch,1,1])),4);
        imf_pha=reshape(imf_pha,[nx,ny,nz,ncontr,nvf]);
        % imf=squeeze(sum(imf_uncomb.*conj(imf_uncomb(:,:,:,:,1,:)),4));
        imf=imf_mag.*exp(1i*angle(imf_pha));
        save_mat(fn_f,'imf',...
                 'overwrite',1,'read4all',1);
        %
        if par.uw
            disp('*** Unwarping full-fov images ***');
            % extract B1
            imf_1=imf_mag(:,:,:,1,:);
            [~,mask_thrd]=mask1d(col(imf_1(:,:,:,1,par.isen)),0.3,0.15);
            b1=sense_m(imf_1,squeeze(imf_1)>mask_thrd,1,'cp_all',1);
            [imf_uw,b0_uw]=ste_uw_inv_mv(imf,data.pe.sp,...
                                         data.te,...
                                         mask_thrd,...
                                         pe_dir,...
                                         1,...
                                         2,...
                                         [],[],b1,1);
            save_mat(fn_f_uw,'imf_uw','b0_uw',...
                     'overwrite',1,'read4all',1);
        end
    end
    % ----------------------------------------------- %
    % reconstruct accelerated volumes for motion estimation
    % ----------------------------------------------- %
    if ismember('motrec',chain)
        % load(fn_cp_ref);
        % tmp=read_data(fn_sen_1v);
        % s_ste_1v=tmp.sensit;
        load(fn_grappa_calib);
        load(fn_sen_1v);
        gcalib.w=w;
        gcalib.par=par_g;
        disp('*** Reconstructing accelerated images ***');
        im_fast=proc_sorted_ste(sorted_ste,'grappa',gcalib);
        im_fast(isnan(im_fast))=0;
        [nx,ny,nz,nch,ncontr,nr,nvf]=size(im_fast);
        im_mot_mag=reshape(sum(abs(im_fast),4),...
                           [nx,ny,nz,ncontr,nr,nvf]);
        im_mot_pha=sum((abs(im_fast(:,:,:,:,1,:,:)).*im_fast).*...
                       exp(-1i*reshape(cp,[1,1,1,nch,1,1,1])),4);
        im_mot_pha=reshape(im_mot_pha,[nx,ny,nz,ncontr,nr,nvf]);
        im_mot=im_mot_mag.*...
               exp(1i*angle(im_mot_pha));
        if field_true(par,'shift_ste')
            im_mot=shift(im_mot,pix_ste_shift);
        end
        save_mat(fn_im_mot,'im_mot','overwrite',1);
        if par.fast_uw
            disp('*** Unwarping fast images ***');
            % extract B1
            im_fast_mag=abs(im_mot(:,:,:,1,:,:));
            im_fast_mag=reshape(im_fast_mag,[nx,ny,nz,1,nr*nvf]);
            [~,mask_thrd]=mask1d(col(im_fast_mag(:,:,:,1,2)),0.3,0.15);
            b1=sense_m(im_fast_mag,squeeze(im_fast_mag)>mask_thrd,1,'cp_all',1);
            im_mot=reshape(im_mot,[nx,ny,nz,ncontr,nr*nvf]);
            [im_mot_uw,b0_mot_uw]=ste_uw_inv_mv(im_mot,data.pe.sp,...
                                                data.te,...
                                                mask_thrd,...
                                                pe_dir,...
                                                1,...
                                                2,...
                                                [],[],b1,1);
            im_mot_uw=reshape(im_mot_uw,[nx,ny,nz,ncontr,nr,nvf]);
            b0_mot_uw=reshape(b0_mot_uw,[nx,ny,nz,nr,nvf]);
            save_mat(fn_im_mot_uw,'im_mot_uw','b0_mot_uw','overwrite',1);
        end
    end
    
    % ----------------------------------------------- %
    % estimate motion using accelerated images
    % ----------------------------------------------- %
    if ismember('motest',chain)
        disp('*** Estimating motion from accelerated images ***');
        load(fn_im_mot);
        load(fn_ste_mask);
        if exist(fn_ste_info)
            load(fn_ste_info);
        else
            nvf_ceil=ceil(sorted_ste.nshot_nblpl/sorted_ste.nshot_fv);
            fcrptd=zeros(nvf_ceil,1);
            ste_info.fcrptd=fcrptd;
        end
        % load(fn_sen_1v);
        nz=size(mask_brain,3);
        [~,~,~,~,~,nv]=size(im_mot);
        if strcmp(vendor,'siemens')
            im_ref=squeeze(abs(im_mot(:,:,:,1,:,par.isen)));
        elseif strcmp(vendor,'philips')
            im_ref=squeeze(abs(im_mot(:,:,:,1,1,par.isen)));
        end
        [mot_par,im_coreg_p]=ste_mot_est(squeeze(im_mot(:,:,:,1,:,:)),...
                                         im_ref,...
                                         mask_brain_sh,...
                                         par.mot_method,...
                                         1,...
                                         res_ste);
        mot_par_raw=mot_par;
        save_mat(fn_motpar_raw,'mot_par_raw',...
                 'overwrite',1,'read4all',1);
        % Certain volumes near the end are not acquired actually
        if nv*sorted_ste.nparv_cyc>sorted_ste.nv
            mot_par(:,sorted_ste.nv+1:nv*sorted_ste.nparv_cyc)=...
                repmat(mot_par(:,sorted_ste.nv),...
                       [1,nv*sorted_ste.nparv_cyc-sorted_ste.nv]);

        end
        % unstable in the first volume
        mot_par(:,:,1)=repmat(mot_par(:,1,2),[1,sorted_ste.nparv_cyc]);
        
        %
        if isfield(par,'regress_nav') && par.regress_nav>0
            
            if isfield(sorted_ste.para,'loop_order') && ...
                    strcmp(sorted_ste.para.loop_order,'zy_order') && ...
                    sorted_ste.para.isgre
                cyc_regress=sorted_ste.para.n_partitions/...
                    sorted_ste.para.sense_rate_s;
            else
                cyc_regress=sorted_ste.para.n_interleaves;
            end
            
            si_mot=size(mot_par);
            nfullnav=si_mot(3);
            nnav=si_mot(2)*si_mot(3);
            mot_par=reshape(mot_par,[6,nnav]).';
            
            for i=1:6
                if fullSTEAcq
                    mot_par_tmp=mot_par(:,i);
                    mot_par_tmp=reshape(mot_par_tmp,...
                                        [sorted_ste.nparv_cyc,nfullnav]);
                    for ir=1:sorted_ste.nparv_cyc
                        mot_par_tmp(ir,:)=regress_harm(mot_par_tmp(ir,:).',...
                                              sorted_ste.para.n_main_tr,...
                                              cyc_regress,...
                                              par.regress_nav,...
                                              sorted_ste.nshot_fv+sorted_ste.nblpl_fv);
                        mot_par_tmp(ir,:)=mot_par_tmp(ir,:)-mot_par_tmp(ir,par.isen);
                        mot_par_tmp(ir,1)=mot_par_tmp(ir,2);
                        mot_par_tmp(ir,end)=mot_par_tmp(ir,end-1);
                    end
                    mot_par(:,i)=mot_par_tmp(:);
                else
                    mot_par(:,i)=regress_harm(mot_par(:,i),...
                                              nnav*(sorted_ste.nshot_v+1),...
                                              cyc_regress,...
                                              par.regress_nav,...
                                              sorted_ste.nshot_v+1);
                end
                
                
            end
            
            mot_par=reshape(mot_par.',si_mot);
        end
        
        % remove the system error due to different fast acquision k-space positions
        if fullSTEAcq
            mot_par=ste_rect_mot_fast(mot_par,par);
        end
        
        save_mat(fn_motpar_nuw,'mot_par',...
                 'overwrite',1,'read4all',1);
        
        % refine determination of corrupted with motion
        fcrptdm=ste_def_crpt_m(mot_par,par,sorted_ste,res_ste);
        

        
        ste_info.fcrptdm=fcrptdm;
        ste_info.fcrptd=fcrptdm | ste_info.fcrptd;
        %
        ste_info.fcrptd(par.isen)=0;
        disp(['*** ' int2str(total(ste_info.fcrptd)) '/' ...
              int2str(length(ste_info.fcrptd)) ...
              ' volumes are corrupted in the end ***']);
        %% unwarpping navigator is seldomly used.
        %% code is commented out for now.
        %% 09042023 by Jiaen Liu
% $$$         if field_true(par,'fast_uw')
% $$$             load(fn_im_mot_uw);
% $$$             if strcmp(vendor,'siemens')
% $$$                 im_ref=squeeze(abs(im_mot_uw(:,:,:,1,:,par.isen)));
% $$$             elseif strcmp(vendor,'philips')
% $$$                 im_ref=squeeze(abs(im_mot_uw(:,:,:,1,1,par.isen)));
% $$$             end
% $$$             [mot_par_uw,im_coreg_p]=ste_mot_est(squeeze(im_mot_uw(:,:,:,1,:,:)),...
% $$$                                                 im_ref,...
% $$$                                                 mask_brain_sh,...
% $$$                                                 par.mot_method,...
% $$$                                                 1,...
% $$$                                                 res_ste);
% $$$             if nv*sorted_ste.nparv_cyc>sorted_ste.nv
% $$$                 mot_par_uw(:,sorted_ste.nv+1:nv*sorted_ste.nparv_cyc)=...
% $$$                     repmat(mot_par_uw(:,sorted_ste.nv),...
% $$$                            [1,nv*sorted_ste.nparv_cyc-sorted_ste.nv]);
% $$$ 
% $$$             end
% $$$             mot_par_uw(:,:,1)=repmat(mot_par_uw(:,1,2),[1,sorted_ste.nparv_cyc]);
% $$$             % remove the system error due to different fast acquision k-space positions
% $$$             if fullSTEAcq
% $$$                 mot_par_uw=ste_rect_mot_fast(mot_par_uw,par);
% $$$             end
% $$$             mot_par=mot_par_uw;
% $$$         end
        
        save_mat(fn_motpar,'mot_par',...
                 'overwrite',1,'read4all',1);
        if field_true(par,'save_coreg')
            save_mat(fn_im_coreg_p,'im_coreg_p',...
                     'overwrite',1,'read4all',1);
        end
        ste_info.mot_par=mot_par;
        save_mat(fn_ste_info,'ste_info',...
                 'overwrite',1,'read4all',1);
    end
    % ----------------------------------------------- %
    % Estimate motion using full-fov images
    % ----------------------------------------------- %
    if ismember('fmotest',chain) && fullSTEAcq
        disp('*** Estimating motion from full-fov images ***');
        load(fn_f);

        load(fn_ste_mask);
        nz=size(mask_brain,3);
        si=size(imf);
        nv=si(5);
        [mot_par_f,im_coreg_f]=ste_mot_est(squeeze(imf(:,:,:,1,:)),...
                                           par.isen,...
                                           mask_brain_sh,...
                                           par.mot_method,...
                                           0,...
                                           res_ste);
        % first volume not stable
        mot_par_f(:,1)=mot_par_f(:,2);
        % the last volume may not be completed.
        mot_par_f(:,end)=mot_par_f(:,sorted_ste.nvf);
        if isfield(par,'regress_nav') && par.regress_nav>0
            
            if isfield(sorted_ste.para,'loop_order') && ...
                    strcmp(sorted_ste.para.loop_order,'zy_order') && ...
                    sorted_ste.para.isgre
                cyc_regress=sorted_ste.para.n_partitions/...
                    sorted_ste.para.sense_rate_s;
            else
                cyc_regress=sorted_ste.para.n_interleaves;
            end
            for i=1:6
                mot_par_f(i,:)=regress_harm(mot_par_f(i,:).',...
                                            sorted_ste.para.n_main_tr,...
                                            cyc_regress,...
                                            par.regress_nav,...
                                            sorted_ste.nshot_fv+sorted_ste.nblpl_fv);
                mot_par_f(i,:)=mot_par_f(i,:)-mot_par_f(i,par.isen);
                mot_par_f(i,1)=mot_par_f(i,2);
                mot_par_f(i,end)=mot_par_f(i,end-1);
            end
        end
        % 
        save_mat(fn_motparf_nuw,'mot_par_f',...
                 'overwrite',1,'read4all',1);
        %% unwarpping navigator is seldomly used.
        %% code is commented out for now.
        %% 09042023 by Jiaen Liu        
% $$$         if par.uw
% $$$             load(fn_f_uw);
% $$$             [mot_par_f,im_coreg_f]=ste_mot_est(squeeze(imf_uw(:,:,:,1,:)),...
% $$$                                                par.isen,...
% $$$                                                mask_brain_sh,...
% $$$                                                par.mot_method,...
% $$$                                                0,...
% $$$                                                res_ste);
% $$$             % first volume not stable
% $$$             mot_par_f(:,1)=mot_par_f(:,2);
% $$$             % the last volume may not be completed.
% $$$             mot_par_f(:,end)=mot_par_f(:,sorted_ste.nvf);
% $$$         end
        % 
        save_mat(fn_motparf,'mot_par_f',...
                 'overwrite',1,'read4all',1);
        if field_true(par,'save_coreg')
            save_mat(fn_im_coreg_f,'im_coreg_f',...
                     'overwrite',1,'read4all',1);
        end
        load(fn_ste_info);
        ste_info.mot_par_f=mot_par_f;
        load(fn_motpar);
        mot_par_comb=...
            ste_comb_motion_beta(ste_info, ...
                                 mot_par_f,...
                                 mot_par,par,sorted_ste);
        
        ste_info.mot_par=mot_par_comb;
        save_mat(fn_motpar_comb,'mot_par_comb',...
                 'overwrite',1,'read4all',1);
        save_mat(fn_ste_info,'ste_info',...
                 'overwrite',1,'read4all',1);
    end

    % ----------------------------------------------- %
    % Calculate B0 field changes in the head frame
    % ----------------------------------------------- %
    if ismember('db0',chain)
        disp('*** Calculating B0 field changes in the head frame ***');
        load(fn_motpar_nuw);
        load(fn_im_mot);
        if fullSTEAcq
            load(fn_f);
            load(fn_motparf_nuw);
        end
        load(fn_ste_info);
        load(fn_ste_mask);
        si=size(im_mot);
        nv=si(end);
        nx=si(1);
        ny=si(2);
        nz=si(3);
        
        db0_p=zeros(nx,ny,nz,sorted_ste.nparv_cyc,nv);
        im_ref_fast=reshape(mean(im_mot(:,:,:,:,:,par.isen),5),...
                            [nx,ny,nz,ncontr_ste]);
        if isfield(par,'mot_par_ref') && ~isempty(par.mot_par_ref)
            mot_par=cat_mot_par(mot_par,par.mot_par_ref);
            mot_par_f=cat_mot_par(mot_par_f,par.mot_par_ref);
        end
        mask_db0_brain=mask_brain;
        if mecho
            dte=(sorted_ste.te(end/2+1)-sorted_ste.te(1))*1e-3;
            % calculate b0 changes using accelerated images
            for i=1:sorted_ste.nparv_cyc
                if strcmp(vendor,'siemens')
                    im_ref=im_mot(:,:,:,:,i,par.isen);
                elseif strcmp(vendor,'philips')
                    im_ref=im_mot(:,:,:,:,1,par.isen);
                end
                print_countdown(sorted_ste.nparv_cyc,...
                                i,...
                                'Number of accelerated types left:');
                mot_par_tmp=squeeze(mot_par(:,i,:));
                if field_true(par,'db0_ignore_mot')
                    mot_par_tmp(:)=0;
                end
                [db0_tmp,im_warp_tmp,mask_in_fov]=ste_db0_headframe(im_mot(:,:,:,:,i,:),...
                                                                  im_ref,...
                                                                  mot_par_tmp,...
                                                                  dte,ncontr_ste,...
                                                                  res_ste/par.ste_intp_res);
                mask_db0_brain=mask_db0_brain&mask_in_fov;
                if field_true(par,'crct_eddy')
                    db0_eddy=ste_db0_eddy(data.para,...
                                          affine_forward(mot_par(:,i,:)),...
                                          par.crct_eddy_fn);
                else
                    db0_eddy=0;
                end
                db0_p(:,:,:,i,:)=db0_tmp+db0_eddy;
            end
            if fullSTEAcq
                % calculate b0 changes using full-fov images
                mot_par_tmp=mot_par_f;
                if field_true(par,'db0_ignore_mot')
                    mot_par_tmp(:)=0;
                end
                [db0_f,imf_warp_tmp,mask_in_fov]=...
                    ste_db0_headframe(imf,...
                                      par.isen,...
                                      mot_par_tmp,...
                                      dte,ncontr_ste,...
                                      res_ste/par.ste_intp_res);
                mask_db0_brain=mask_db0_brain&mask_in_fov;
                if field_true(par,'crct_eddy')
                    db0_eddy=ste_db0_eddy(data.para,...
                                          affine_forward(mot_par_f),...
                                          par.crct_eddy_fn);
                else
                    db0_eddy=0;
                end
                db0_f=db0_f+db0_eddy;
            end
% $$$             end
        else
            necho_contr=sorted_ste.necho/sorted_ste.ncontr;
            dte=sorted_ste.te(floor(necho_contr/2)+1)*1e-3;
            im_mot_1=im_mot(:,:,:,1,:,:);
            
            for i=1:sorted_ste.nparv_cyc
                if strcmp(vendor,'siemens')
                    im_ref=im_mot_1(:,:,:,:,i,par.isen);
                elseif strcmp(vendor,'philips')
                    im_ref=im_mot_1(:,:,:,:,1,par.isen);
                end
                print_countdown(sorted_ste.nparv_cyc,...
                                i,...
                                'Number of accelerated types left:');
% $$$                 tra_matrix_tmp=motpar2affine(squeeze(mot_par(:,i,:)),...
% $$$                                              [nx,ny,nz],'inverse',...
% $$$                                              res_ste/par.ste_intp_res);
                mot_par_tmp=squeeze(mot_par(:,i,:));
                if field_true(par,'db0_ignore_mot')
                    mot_par_tmp(:)=0;
                end
                [db0_tmp,im_warp_tmp,mask_in_fov]=...
                    ste_db0_headframe(im_mot_1(:,:,:,:,i,:),...
                                      im_ref,...
                                      mot_par_tmp,...
                                      dte,1,...
                                      res_ste/par.ste_intp_res);
                mask_db0_brain=mask_db0_brain&mask_in_fov;
                if field_true(par,'crct_eddy')
                    db0_eddy=ste_db0_eddy(data.para,...
                                          affine_forward(mot_par(:,i,:)),...
                                          rp('mid1702.db0_ste_gre_coef.mat','20180808_1'));
                else
                    db0_eddy=0;
                end
                db0_p(:,:,:,i,:)=db0_tmp+db0_eddy;
            end
            if fullSTEAcq
                mot_par_tmp=mot_par_f;
                if field_true(par,'db0_ignore_mot')
                    mot_par_tmp(:)=0;
                end
                [db0_f,imf_warp_tmp,mask_in_fov]=...
                    ste_db0_headframe(imf(:,:,:,1,:),...
                                      par.isen,...
                                      mot_par_tmp,...
                                      dte,1,...
                                      res_ste/par.ste_intp_res);
                mask_db0_brain=mask_db0_brain&mask_in_fov;
                if field_true(par,'crct_eddy')
                    db0_eddy=ste_db0_eddy(data.para,...
                                          affine_forward(mot_par_f),...
                                          rp('mid1702.db0_ste_gre_coef.mat','20180808_1'));
                else
                    db0_eddy=0;
                end
                db0_f=db0_f+db0_eddy;
            end
        end
        % regress out phase encoding of the main sequence related artifact
        if isfield(par,'regress_nav') && par.regress_nav>0
            
            if isfield(sorted_ste.para,'loop_order') && ...
                    strcmp(sorted_ste.para.loop_order,'zy_order') && ...
                    sorted_ste.para.isgre
                cyc_regress=sorted_ste.para.n_partitions/...
                    sorted_ste.para.sense_rate_s;
            else
                cyc_regress=sorted_ste.para.n_interleaves;
            end
            if fullSTEAcq
                % fast navigator
                for ir=1:sorted_ste.nparv_cyc
                    db0_p_tmp=squeeze(db0_p(:,:,:,ir,:));
                    si_db0=size(db0_p_tmp);
                    db0_p_tmp=combine_dim(db0_p_tmp,[1,2,3]).';
                    db0_p_tmp=regress_harm(db0_p_tmp,...
                                           sorted_ste.para.n_main_tr,...
                                           cyc_regress,...
                                           par.regress_nav,...
                                           sorted_ste.nshot_fv+sorted_ste.nblpl_fv);
                    db0_p_tmp=db0_p_tmp-db0_p_tmp(par.isen,:);
                    db0_p_tmp(1,:)=db0_p_tmp(2,:);
                    db0_p_tmp(end,:)=db0_p_tmp(end-1,:);
                    db0_p_tmp=reshape(db0_p_tmp.',si_db0);
                    db0_p(:,:,:,ir,:)=db0_p_tmp;
                end
                % full navigator
                si_db0_f=size(db0_f);
                db0_f_tmp=combine_dim(db0_f,[1,2,3]).';
                db0_f_tmp=regress_harm(db0_f_tmp,...
                                       sorted_ste.para.n_main_tr,...
                                       cyc_regress,...
                                       par.regress_nav,...
                                       sorted_ste.nshot_fv+sorted_ste.nblpl_fv);
                db0_f_tmp=db0_f_tmp-db0_f_tmp(par.isen,:);
                db0_f_tmp(1,:)=db0_f_tmp(2,:);
                db0_f_tmp(end,:)=db0_f_tmp(end-1,:);
                db0_f=reshape(db0_f_tmp.',si_db0_f);
            else
                si_db0=size(db0_p);
                nnav=si_db0(4)*si_db0(5);
                db0_p_tmp=combine_dim(combine_dim(db0_p,[4,5]),[1,2,3]).';
                db0_p_tmp=regress_harm(db0_p_tmp,...
                                           nnav*(sorted_ste.nshot_v+1),...
                                           cyc_regress,...
                                           par.regress_nav,...
                                           sorted_ste.nshot_v+1);
                db0_p=reshape(db0_p_tmp.',si_db0);
            end
        end
        
        % rectify db0_p because of effect from fast acquisition
        if fullSTEAcq
            db0_p=...
                ste_rect_b0_fast(db0_p,par);
        end
        
        if fullSTEAcq
            save_mat(fn_db0_f,'db0_f',...
                     'overwrite',1,'read4all',1);
        end
        if field_true(par,'save_db0')
            save_mat(fn_db0_p,'db0_p',...
                 'overwrite',1,'read4all',1);
        end
        load(fn_ste_info);
        ste_info.mask_db0_brain=mask_db0_brain;
        save_mat(fn_ste_info,'ste_info',...
                 'overwrite',1,'read4all',1);
        % combine fast and full navigators
        if fullSTEAcq
            db0_comb=ste_comb_b0(...
                db0_f,db0_p,ste_info,par,sorted_ste);
        else
            db0_comb=db0_p;
        end
        % temporal and spatial filtering
        if isfield(par,'db0_comb_filter') && ~isempty(par.db0_comb_filter)
            db0_comb_filter=par.db0_comb_filter;
            si=size(db0_comb);
            db0_comb=combine_dim(db0_comb,[4,5]);
            h=ndgauss(round_odd(db0_comb_filter*6),db0_comb_filter);
            db0_comb=imfilter(db0_comb,h);
            db0_comb=reshape(db0_comb,si);
        end
        save_mat(fn_db0_comb,'db0_comb',...
                 'overwrite',1,'read4all',1);
    end
    
    % ----------------------------------------------- %
    % cluster data based on smoothed motion
    % for each cluster a sensivitiy map
    % will be created later
    % ----------------------------------------------- %
    if ismember('cluster',chain)
        load(fn_ste_info);
        disp('*** Clustering the GRE data into groups ***');
        if field_true(par,'use_pri_mot')
            load(file_addext(fn_motpar,'_pri'));
        else
            load(fn_motpar);
        end
% $$$         if fullSTEAcq
% $$$             load(fn_motparf);
% $$$         end
        load(fn_ste_mask);
        load(fn_db0_comb);
        if isfield(par,'mot_par_ref') && ~isempty(par.mot_par_ref)
            mot_par=cat_mot_par(mot_par,par.mot_par_ref);
            mot_par_f=cat_mot_par(mot_par_f,par.mot_par_ref);
        end
        para_tmp=sorted_ste.para;
        para_tmp.idx_shot_cenpf=sorted_ste.idx_shot_cenpf;
        para_tmp.idx_shot_cenff=sorted_ste.idx_shot_cenff;
        para_tmp.ste_intp_res=par.ste_intp_res;
        [ste_info,b0_fit]=ste_cluster(ste_info,par,...
                                      para_tmp,...
                                      sorted_ste,db0_comb,ste_info.mask_db0_brain);
        par.idx=ste_info.idx;
        par.nc=ste_info.nc;
        if ~isempty(b0_fit)
            save_mat(fn_b0_fit_cluster,'b0_fit',...
                     'overwrite',1,'read4all',1);
        end
        disp(['*** There are ' int2str(ste_info.nc) ' clusters. ***']);
        save_mat(fn_ste_info,'ste_info',...
                 'overwrite',1,'read4all',1);
    end


    % ----------------------------------------------- %
    % Final reconstruction
    % ----------------------------------------------- %
    % back up on biowulf because the reconstrucion may take too long
    % This is not used any more because the data result folder is not fixed.
    % backup_biowulf_lcl();
    if ismember('grerec',chain)
        disp('*** Preparing the parameters for reconstruction ***');
        % only need to load motparf
        % motpar is stored in ste_info
        if fullSTEAcq
            load(fn_motparf);
            if isfield(par,'mot_par_ref') && ~isempty(par.mot_par_ref)
                mot_par_f=cat_mot_par(mot_par_f,par.mot_par_ref);
            end
        else
            mot_par_f=[];
        end
        load(fn_db0_comb);
        load(fn_ste_mask);
        load(fn_ste_info);
        % convert the affine matrix to matlab version
        nx=data.para.steref_dim_r;
        ny=data.para.steref_dim_p;
        nz=data.para.steref_dim_s;

        ste_info=ste_motion(ste_info,sorted_ste,par);
        % pick uncombined full-fov data for sensivity maps
        if fullSTEAcq
            isen=zeros(ste_info.nc,1);
            for i=1:ste_info.nc
                isen(i)=ste_info.md{i}.isen;
            end
            imf_uncomb=get_k_ste(sorted_ste,isen,nshot_fv);
            sitmp=size(imf_uncomb);
            imf_uncomb=fftmr(imf_uncomb,-1,[1,2,3])*prod(sitmp(1:3));
        elseif strcmp(vendor,'philips') && ...
                isempty(par.use_ext_ref)
            isen=1;
            imf_uncomb=get_k_ste(sorted_ste,isen,nshot_fv);
            sitmp=size(imf_uncomb);
            imf_uncomb=fftmr(imf_uncomb,-1,[1,2,3])*prod(sitmp(1:3));
            mot_par_f=zeros(6,1);
        else
            imf_uncomb=[];
        end
        [para,ste_info,b0_fit,mot_par_test,b1]=...
            ste_final_prep_nintp(sorted_ste,ste_info,mask,...
                                 ste_info.mask_db0_brain,...
                                 db0_comb,...
                                 mot_par_f,imf_uncomb,par);
        save_mat(fn_b0_fit,'b0_fit',...
                 'overwrite',1,'read4all',1);
        save_mat(fn_b1,'b1',...
                 'overwrite',1,'read4all',1);
        % no b0 correction
        if par.nb0
            para.b0(:)=0;
            for i=1:para.nm
                para.md(i).gb0(:)=0;
                para.md(i).db0(:)=0;
            end
        end
        if field_true(par,'n_gb0')
            for i=1:para.nm
                para.md(i).gb0(:)=0;
            end
        end 
        if field_true(par,'n_gb0_db0')
            for i=1:para.nm
                para.md(i).gb0(:)=0;
                para.md(i).db0(:)=0;
            end
        end
        if par.ncorrection
            para.b0(:)=0;
            for i=1:para.nm
                para.md(i).m(:)=0;
                para.md(i).gb0(:)=0;
                para.md(i).db0(:)=0;
            end
        end
        if field_true(par,'nmoco')
            for i=1:para.nm
                para.md(i).m(:)=0;
            end
        end
% $$$         if field_true(par,'use_pri_b0')
% $$$             fn_gre=fullfile(subdir,['mid' int2str(mid) '.k_nav0_pri_uncomb.svd']);
% $$$         end
        nk_shot=para.nk_shot;
        n_partitions=para.n_partitions;
        s1=para.sense_rate_p;
        s2=para.sense_rate_s;
        nshot=para.np*n_partitions/...
              s1/s2/nk_shot;
        idx_reorder_shot=zeros(nshot,1);
        j=0;
        if par.discard
            % only use data from stable full-fov ste acqs
            % create motion mask for the option of discarding data
            for i=1:para.nm
                n_shot_cluster=para.md(i).nk/nk_shot;
                if mod(n_shot_cluster,1)~=0
                    error('*** Number of shots per cluster is not an integer! ***');
                end
                idx_reorder_shot(j+1:j+n_shot_cluster)=...
                    find(ste_info.idx_gre==i&...
                         ste_info.mask_gre);
                j=j+n_shot_cluster;
            end
        else
            for i=1:para.nm
                n_shot_cluster=para.md(i).nk/nk_shot;
                if mod(n_shot_cluster,1)~=0
                    error('*** Number of shots per cluster is not an integer! ***');
                end
                idx_reorder_shot(j+1:j+n_shot_cluster)=find(ste_info.idx_gre==i);
                j=j+n_shot_cluster;
            end
        end
        idx_reorder_shot(j+1:end)=[];
        para.idx_reorder_shot=idx_reorder_shot;
        % filter B0
        h=fgaussian3([7,7,7],[0.6,0.6,0.6]);
        for i=1:ste_info.nc
            para.b0(:,:,:,i)=imfilter(para.b0(:,:,:,i),h);
        end
        db0_c=para.b0;
        save_mat(fn_db0_cluster,'db0_c',...
                 'overwrite',1,'read4all',1);
        save_mat(fn_ste_info,'ste_info',...
                 'overwrite',1,'read4all',1);
        disp('*** Reconstructing GRE images ***');
        %
        ncontr=length(par.icontr);
        kd=double(read_data(fn_gre));

        if kd==-1
            error('*** K space data does not exist! ***');
        end
        if strcmp(vendor,'siemens')
            kd=reshape(kd,[para.nr,para.np/para.sense_rate_p,...
                           para.n_partitions/para.sense_rate_s,...
                           para.n_slices,nch,...
                           necho,para.n_reps]);
            kd=kd(:,:,:,:,:,par.icontr,:);
            if strcmp(par.precision,'single')
                para.b1n=single(para.b1n);
                kd=single(kd);
            end
            if ~para.isgre
                % if reconstructing EPI, need to swap the 
                % EPI dimension
                kd=reshape(kd,[para.nr,n_interl,nk_shot,npar_acc,...
                               para.n_slices,nch,ncontr,...
                               para.n_reps]);
                kd=permute(kd,[1,3,2,4,5,6,7,8]);
            end
        elseif strcmp(vendor,'philips')
            kd=reshape(kd,[para.nr,para.nk_shot,...
                           para.nte_contr,para.n_slices,nch,...
                           n_interl*npar_acc,para.n_reps]);
            kd=kd(:,:,par.icontr,:,:,:,:);
            kd=permute(kd,[1,2,6,4,5,3,7]);
        end
        kd=reshape(kd,[para.nr,nk_shot,n_interl*npar_acc,...
                       para.n_slices,nch,ncontr,nreps]);
        if par.simul_noise
            % divide cov by 2 because of oversampling
            sqr_cov=chol(inv(sorted_ste.inv_cov)/2)';
            kd=sqr_cov*randn(nch,para.nr*nk_shot*n_interl*npar_acc*...
                                para.n_slices*ncontr*nreps);
            kd=reshape(kd,[nch,para.nr,nk_shot,n_interl*npar_acc,...
                           para.n_slices,ncontr,nreps]);
            kd=permute(kd,[2,3,4,5,1,6,7]);
        end
        % rearrange the k-space order to match 
        % the order defined in the clustering
        % algorithm
        if ste_info.nc>1
            kd=kd(:,:,idx_reorder_shot,:,:,:,:);
        end
        % start to drop slice as a dimention because this code is meant for 
        % 3D acqusition with one slab(slice)
        kd=reshape(kd,[numel(kd)/nch/ncontr/nreps,nch,ncontr,nreps]);
        sqr_inv_cov=conj(chol(sorted_ste.inv_cov,'lower'));
        for irep=1:nreps
            for ic=1:ncontr
                kd(:,:,ic,irep)=kd(:,:,ic,irep)*sqr_inv_cov;
            end
        end
        %% down sample a full fov dataset
        if isfield(par,'dsamp')
            % check to make sure there is no discard, no acceleration
            if field_true(par,'discard') || ...
                    para.sense_rate_p>1 || ...
                    para.sense_rate_s>1 || ...
                    ~para.isgre
                warning('*** Downsampling was not performed because it''s not fully sampled or the discard option is on or it''s EPI! ***');
            else
                [para,kd,par]=ste_dsamp(para,kd,par);
            end
        end
        para.max_nch_seg=par.max_nch_seg;
        para.pcaflag=par.pcaflag;
        para.npca=par.npca;
        para.n_seg_mem=par.n_seg_mem;
        if field_true(par,'pcaflag')
            pca_struct.coefsinterp=ste_info.coefsinterp;
            pca_struct.scorefilt=ste_info.scorefilt;
        end
        clearvars -except para_mr para kd par ...
            ncontr nreps pca_struct;
        recon_t=zeros(ncontr,nreps);
        epsil=[];
        im_recon=zeros(para.nr,para.np,...
                       para.n_partitions*para.n_slices,...
                       ncontr,nreps,par.precision);

        nshot_rep=para.n_interleaves*para.n_partitions/para.sense_rate_s;
        nk_shot=para.nk_shot;
        for irep=1:nreps
            disp(['*** Reconstructing repetition ' ...
                  int2str(irep) ...
                  ' of ' int2str(nreps) ...
                  ' ***']);
            idx_shot_rep=[1:nshot_rep]+(irep-1)*nshot_rep;
            for ic=1:ncontr

                disp(['*** Reconstructing contrasts ' ...
                      int2str(ic) ...
                      ' of ' int2str(length(par.icontr)) ...
                      ' ***']);
                if par.pcaflag
                    % retrieve coeffs for the current repetition
                    para.coefsinterp=...
                        pca_struct.coefsinterp(idx_shot_rep,:,par.icontr(ic));
                    % duplicate the coeffs for the k-space lines in EPI case
                    para.coefsinterp=reshape(para.coefsinterp,...
                                             [1,nshot_rep,par.npca]);
                    para.coefsinterp=repmat(para.coefsinterp,...
                                            [nk_shot,1,1]);
                    para.coefsinterp=reshape(para.coefsinterp,...
                                             [nk_shot*nshot_rep,par.npca]);
                    para.scorefilt=...
                        pca_struct.scorefilt(:,:,par.icontr(ic));
                end
                [im_recon(:,:,:,ic,irep),flag,epsil_cur,b,~,~,t]=...
                    recon_mb_epi_beta(kd(:,:,ic,irep),...
                                      para,[par.icontr(ic),irep],...
                                      par.precision);
                recon_t(ic,irep)=t;
                epsil=[epsil;epsil_cur(:)];
                im_tmp=im_recon(:,:,:,ic,irep);
                if field_true(par,'save_temp')
                    save im_recon_tmp.mat im_tmp ic irep;
                    siem_to_nifti('im_recon_tmp.nii.gz',abs(im_tmp),para_mr,1,1);
                end
            end
        end
        par.recon_t=recon_t;
        par.hostname=getenv('HOSTNAME');
        par.epsil=epsil;
        par.flag=flag;
        par.para=para_mr;
    end
end

