% version
% v1.1: add apodization on 08/19/2019
% Translated from IDL to Matlab on 09/25/2020
%
% 2022-09-29 PvG & JAdZ: Added -filemode 0666 to sort_siemens call to allow all to overwrite files
% 2022-10-05 Jiaen Liu: Fix linear phase issue in recon
% 2022-12-20 Jiaen Liu: Philips support
% 2023-05 Jiaen Liu: Philips sense reference support
% 2024-08-18 Jiaen Liu: address channel mismatch between reference and main scans
function s=prep_ste(mid,varargin)
    version = 'v1.1';
    p=inputParser;
    nstd_pc=1;
    frac_pc=0.3;
    no_pc=0;
    no_b0_main=0;
    apodization=0.3;
    no_save=0;
    no_main=0;
    mid_pimg=[];
    mid_blpo=[];
    data_path='.';
    vendor='siemens';
    sense_philips=0;
    regr_order=0;
    pe_path='/data/user/jiaen_liu/pe_files/';
    if contains(getenv('HOSTNAME'),'bobo','IgnoreCase',true)
        pe_path='/raid/common/fmrif7t/pe_files/';
    elseif contains(getenv('HOSTNAME'),'hawaii','IgnoreCase',true)
        pe_path='/data/user/jiaen_liu/pe_files/';
    end
    addParameter(p,'nstd_pc',nstd_pc,@isnumeric);
    addParameter(p,'frac_pc',frac_pc,@isnumeric);
    addParameter(p,'no_pc',no_pc,@isnumeric);
    addParameter(p,'no_b0_main',no_pc,@isnumeric);
    addParameter(p,'apodization',apodization,@isnumeric);
    addParameter(p,'no_save',no_save,@isnumeric);
    addParameter(p,'no_main',no_main,@isnumeric);
    addParameter(p,'mid_pimg',mid_pimg,@(x)isnumeric(x)||ischar(x));
    addParameter(p,'mid_blpo',mid_blpo,@(x)isnumeric(x)||ischar(x));
    addParameter(p,'data_path',data_path,@ischar);
    addParameter(p,'pe_path',pe_path,@ischar);
    addParameter(p,'vendor',vendor,@ischar);
    addParameter(p,'sense_philips',sense_philips,@isnumeric);
    addParameter(p,'regr_order',regr_order,@isnumeric);
    p.parse(varargin{:});
    nstd_pc=p.Results.nstd_pc;
    frac_pc=p.Results.frac_pc;
    no_pc=p.Results.no_pc;
    no_b0_main=p.Results.no_b0_main;
    apodization=p.Results.apodization;
    no_save=p.Results.no_save;
    no_main=p.Results.no_main;
    mid_pimg=p.Results.mid_pimg;
    mid_blpo=p.Results.mid_blpo;
    data_path=p.Results.data_path;
    pe_path=p.Results.pe_path;
    vendor=p.Results.vendor;
    sense_philips=p.Results.sense_philips;
    regr_order=p.Results.regr_order;
    if numel(mid_pimg)>1
        % 2024-08-18
        error('*** Up to one refererence scan is supported! ***');
    end
    % mid is a char in standard alone application
    if ischar(mid)
        mid=eval(mid);
    end
    if ischar(mid_pimg)
        mid_pimg=eval(mid_pimg);
    end
    if ischar(mid_blpo)
        mid_blpo=eval(mid_blpo);
    end
    n = length(mid);
    cd(data_path);
    
    for imid=1:n
        fn_out=rp(['mid',num2str(mid(imid)),'.steref4recon.svd']);
        if strcmp(vendor,'siemens')
            disp(['*** MID:', num2str(mid(imid)), ' ***']);
            disp('*** Processing navigator data ***');
            para=extract_para(mid(imid));
            fname=get_file_filter('.',['MID*',num2str(mid(imid)),'.nav.svd']); % used to be steref, PvG 21Sep22
            if isempty(fname)
                fname=get_file_filter('.',['MID*',num2str(mid(imid)),'.steref.svd']); % legacy support for older data
                if isempty(fname)
                    cmd=['sort_siemens -filemode 0666 ' num2str(mid(imid))];
                    if system(cmd)~=0
                        error(['*** ',cmd,' was not successful! ***']);
                    end
                    fname=get_file_filter('.',['MID*',num2str(mid(imid)),'.nav.svd']); % new sort_siemens output
                    if isempty(fname)
                        fname=get_file_filter('.',['MID*',num2str(mid(imid)),'.steref.svd']); % legacy support for older data
                    end
                end
            end
            fnmdh=get_file_filter('.',['MID*',num2str(mid(imid)),'.mdh']);
        end
        if strcmp(vendor,'siemens')
            d=read_data(fname);
            mask_channel=ones(para.n_channels,1);
            if ~isempty(mid_pimg)
                para_b1=extract_para(mid_pimg);
                idx_channel=[1:para.n_channels].';
            else
                idx_channel=[];
            end
            
        elseif strcmp(vendor,'philips')
            [d,pe,im_nav]=recon_bmir_nav3d(mid(imid));
            para=extract_para_philips(mid(imid));
            if ~isempty(mid_pimg)
                para_b1=extract_para_philips(mid_pimg);
                [idx_channel,mask_channel]=idx_channel_philips(para_b1,para);
            else
                idx_channel=[];
                mask_channel=ones(para.n_channels,1);
            end
        end
        si=size(d);
        necho_ste=si(2);
        nshot=si(end);
        ishot_ref=floor(nshot/2);
        nr=para.steref_dim_r;
        nr_os=nr*2;
        nch=si(end-1);
        nnav=0;
        if para.b_nav_en
            nnav = conditional(para.nav_type==0,1,para.nav_type);
        end
        necho_main=0;
        for icontr=1:para.n_contrasts
            necho_main=necho_main+...
                (para.np/para.n_interleaves/para.sense_rate_p+para.n_refs(icontr));
        end
        % total echos in one shot
        nechot = nnav+necho_ste+necho_main;
        if strcmp(vendor,'siemens')
            % read mdh;
            if contains(para.idea_v,'VB','IgnoreCase',true)
                mdh=defMDH17();
                nByteMdh=size_of(mdh);
                dmdh=read_raw(fnmdh,[nByteMdh,nch,nechot,nshot],'uint8',1);
            elseif contains(para.idea_v,'VD','IgnoreCase',true) || ...
                    contains(para.idea_v,'VE','IgnoreCase',true)
                mdh=defMDH11();
                nByteMdh=size_of(mdh);
                dmdh=read_raw(fnmdh,[nByteMdh,1,nechot,nshot],'uint8',1);
            end
            % STE pe encoding definition
            if para.ste3d_mode==1 || para.ste3d_mode==0 || para.ste3d_mode==4
                % defined by the sequence
                if para.ste3d_mode==1 || para.ste3d_mode==4
                    % zig zag
                    nkyz= para.steref_dim_p*(para.steref_dim_s+1);    
                elseif para.ste3d_mode==0
                    % linear
                    error('*** Linear navigator encoding not supported by the reconstruction! ***');
                    nkyz= para.steref_dim_p*para.steref_dim_s;
                    
                end
                kyz=zeros(2,nkyz);
                if contains(para.idea_v,'VB','IgnoreCase',true)
                    % hard coded zigzap pattern for VB
                    if para.ste3d_mode==1
                        ds=1;
                    elseif para.ste3d_mode==4
                        ds=0;
                    end
                    kyz=[zeros(2,necho_ste),...
                         gen_pe_sense(para.steref_dim_p,para.steref_dim_s,...
                                      4,2,ds,0,0),...
                         gen_pe_sense(para.steref_dim_p,para.steref_dim_s,...
                                      4,2,ds,2,0),...
                         zeros(2,necho_ste),...
                         gen_pe_sense(para.steref_dim_p,para.steref_dim_s,...
                                      4,2,ds,1,0),...
                         gen_pe_sense(para.steref_dim_p,para.steref_dim_s,...
                                      4,2,ds,3,0),...
                         zeros(2,necho_ste),...
                         gen_pe_sense(para.steref_dim_p,para.steref_dim_s,...
                                      4,2,ds,0,1),...
                         gen_pe_sense(para.steref_dim_p,para.steref_dim_s,...
                                      4,2,ds,2,1),...
                         zeros(2,necho_ste),...
                         gen_pe_sense(para.steref_dim_p,para.steref_dim_s,...
                                      4,2,ds,1,1),...
                         gen_pe_sense(para.steref_dim_p,para.steref_dim_s,...
                                      4,2,ds,3,1)];
                else
                    for itr=1:nkyz/necho_ste      
                        for ie=1:necho_ste
                            mdhtmp=cast2struct(dmdh(:,1,nnav+ie,itr+...
                                                    para.n_noise_tr+...
                                                    para.n_blipoff_tr),mdh);
                            kyz(1,ie+ (itr-1)*necho_ste)= mdhtmp.ky;
                            kyz(2,ie+ (itr-1)*necho_ste)= mdhtmp.kz;
                        end
                    end
                    if min(kyz(1,:))>=0
                        kyz(1,:)=kyz(1,:)-floor(para.steref_dim_p/2);
                    end
                    if min(kyz(2,:))>=0
                        kyz(2,:)=kyz(2,:)-floor(para.steref_dim_s/2);
                    end
                end
                pe= struct('hf',zeros(6,1,'single'),...
                           'hl',zeros(2,1,'int32'),...
                           'k',kyz);
                pe.hl(1)=nkyz;
                if para.ste3d_mode==1 || para.ste3d_mode==4
                    pe.y_cyc=2;
                    pe.z_cyc=1;
                else
                    pe.y_cyc=0;
                    pe.z_cyc=0;
                end
                pe.r=8;
                pe.sp=4;
                if para.ste3d_mode==1
                    pe.dkz=1;
                elseif para.ste3d_mode==4
                    pe.dkz=0;
                end
                pe.blipless=2;
                pe.ncontr=necho_ste/para.n_echo_steref;    
            else
                % pe file
                pe_str=typecast(int32(para.ste3d_mode),'uint8');
                pe_str=char(pe_str(1:2));
                pe_file=get_file_filter(pe_path,['pe3d_',pe_str,'_*',...
                                    num2str(necho_ste),'.flt']);
                if ~strcmp(class(pe_file),'char') && numel(pe_file)~=1
                    error('*** The pe_file name is not correct! ***');
                end
                pe=read_pe_files(fullfile(pe_path,pe_file));
            end
        elseif strcmp(vendor,'philips')
            % pe define
            pe_strt.y_cyc=0;
            pe_strt.z_cyc=0;
            pe_strt.blipless=2;
            pe_strt.sp=para.steref_s1;
            r=para.steref_s1*para.steref_s2;
            pe_strt.r=r;
            pe_strt.dkz=floor(para.steref_s2/2);
            pe_strt.ncontr=necho_ste/para.n_echo_steref;
            n_shot_fast_nav3d=para.steref_dim_s/para.steref_s2+1;
            n_shot_full_nav3d=n_shot_fast_nav3d*r;
            pe_strt.k=pe(:,1:n_shot_full_nav3d*necho_ste);
            pe_strt.hl=int32([n_shot_full_nav3d*necho_ste,0]);
            pe=pe_strt;
        end
        % readout direction and te of short reference
        ro_pol=zeros(necho_ste,1);
        te=zeros(necho_ste,1);
        if strcmp(vendor,'siemens')
            for ie=1:necho_ste
                mdhtmp=cast2struct(dmdh(:,1,nnav+ie,1),mdh);
                % refer to !SIEMENS_EVALINFOMASK
                % DEFINE_EVALINFOMASK in EVALINFOMASK_BITS.pro
                ro_pol(ie)=evalmaskbit(mdhtmp,25);
                te(ie)=mdhtmp.te;
            end
            % reverse the data for negative readouts
            d(:,find(ro_pol==0),:,:,:)=flipdim(d(:,find(ro_pol==0),:,:,:),1);
        else
            ro_pol(1:2:end)=1;
            te=para.te_ste(:);
        end
        % regridding
        % for Philips, there is no ramp sampling
        if para.ste_rsamp> 0.01
            % ramp sampling
            mdhtmp=cast2struct(dmdh(:,1,nnav+1,1),mdh);
            para_rsamp=struct('nr',para.steref_dim_r,...
                              't_dwell',para.t_dwell_ste,...
                              'ramp_dur',mdhtmp.ramp,...
                              'ramp_samp_frac',para.ste_rsamp,...
                              'nr_os',para.steref_dim_r*2.0);
            d=regridding_arr(d,para_rsamp,apodization);
        else
            if apodization~=-1
                d=apodize_arr(d,apodization,1);
            end
            d=fftmr(d,-1,1);
        end
        if strcmp(vendor,'siemens')
            d=d(idx_truncate(nr_os,nr),:,:,:,:);
        end
        % remove noise scan
        if strcmp(vendor,'siemens')
            if para.n_noise_shots>0
                d=d(:,:,:,:,para.n_noise_shots+1:end);
                nshot=nshot-para.n_noise_shots;
            end
            % remove blipoff scans
            if para.n_blipoff_reps>0
                n_shot_blipoff=para.n_blipoff_tr;
                d=d(:,:,:,:,n_shot_blipoff+1:end);
                nshot=nshot-n_shot_blipoff;
            else
                n_shot_blipoff=0;
            end
        end
        if strcmp(vendor,'siemens')
            d=reshape(d(:,:,1,:,:),[nr,necho_ste, nch, nshot]);
        end
        if strcmp(vendor,'siemens')
            % odd-even phase different correction
            if ~no_pc
                d=pha_crct_ste(d,pe,nstd_pc,frac_pc);
            end
        end
        % transform to k-space
        d=fftmr(d,1,1);
        % correct TE shift-caused eddy current
        % experimental, not working so far
% $$$       if para.int_te_shift
% $$$           d=steTESClean(d,double(para.n_interleaves),...
% $$$               reshape(pe.k,[2,pe.hl(1)]),para.tr);
% $$$       end
% correct fov shift
        if strcmp(vendor,'siemens')
            d=fov_crct_ste(d,pe.k,para);
        end
        if strcmp(vendor,'siemens')
            % calculate covariance matrix
            cov_mat=covSiem(mid(imid));
        elseif strcmp(vendor,'philips')
            if sense_philips
                % this is rare
                % para_b1.channel_id=para_b1.channel_id(3:end);
                para_b1.channel_id=para_b1.channel_id(idx_channel);
                para_b1.n_channels=length(para_b1.channel_id);
                para_b1.n_slices=1;
            end
% $$$             idx_channel=zeros(length(para.channel_id),1);
% $$$             for i=1:length(para.channel_id)
% $$$                 tmp=find(para_b1.channel_id==...
% $$$                          para.channel_id(i));
% $$$                 if isempty(tmp)
% $$$                     idx_channel(i)=-1;
% $$$                 else
% $$$                     idx_channel(i)=tmp;
% $$$                 end
% $$$             end
            n=read_raw_philips(mid(imid),'type',5);
            % n=n(:,find(idx_channel~=-1));
            n=n(:,find(mask_channel));
            para.n_channels=total(mask_channel);
            cov_mat=cov(conj(n));
        end
        if strcmp(vendor,'siemens')
            if ~contains(para.idea_v,'VB','IgnoreCase',true)
                % get ky and kz for main acquisition
                kyz=zeros(2,necho_main,nshot);
                for i=1:nshot
                    for j=1:necho_main
                        mdhtmp=cast2struct(dmdh(:,1,nnav+necho_ste+j,...
                                                i+n_shot_blipoff+para.n_noise_shots),...
                                           mdh);
                        kyz(1,j,i)=mdhtmp.ky;
                        kyz(2,j,i)=mdhtmp.kz;
                    end
                end
                if min(col(kyz(1,:,:)))>=0
                    kyz(1,:,:)=kyz(1,:,:)-floor(para.np/2);
                end
                if min(col(kyz(2,:,:)))>=0
                    kyz(2,:,:)=kyz(2,:,:)-floor(para.n_partitions/2);
                end
            else
                kyz=[];
            end
        elseif strcmp(vendor,'philips')
            kyz=get_pe_philips(mid(imid),'type',1,'mix',0);
        end
        %
        para.n_channels=total(mask_channel);
        s=struct('d',d(:,:,find(mask_channel),:),...
                 'pe',pe,...
                 'ro_pol',ro_pol,...
                 'te',te,...
                 'para',para,...
                 'cov_mat',cov_mat,...
                 'ishot_ref',ishot_ref,...
                 'kyz',kyz,...
                 'version',version,...
                 'vendor',vendor);
        % remove average B0 in STE
        if strcmp(vendor,'philips')
            regr_order=20;
        end
        if regr_order==0
            par_regress.enable=0;
            par_regress.order=0;
        else
            par_regress.enable=1;
            par_regress.order=regr_order;
        end
        % remove global B0 fluctuation in navigator
        s=detrend_ste(s,ishot_ref,par_regress);
        % process parallel imaging reference data
        if ~isempty(mid_pimg)
            for imid_pimg=1:length(mid_pimg)
                if strcmp(vendor,'siemens')
                    % para_pimg_tmp=extract_para(mid_pimg(imid_pimg));
                    k_pimg_tmp=squeeze(recon_amri_epi(mid_pimg(imid_pimg),'k_return',1));
                elseif strcmp(vendor,'philips')
                    % para_pimg_tmp=extract_para_philips(mid_pimg(imid_pimg));
                    if ~sense_philips
                        k_pimg_tmp=squeeze(recon_bmir_epi(mid_pimg(imid_pimg)));
                    else
                        % this is rare
                        k_pimg_tmp=sense_ref_philips(mid_pimg(imid_pimg));
                    end
                    % nx x np x nslice x nch x necho
                    k_pimg_tmp=k_pimg_tmp(:,:,:,idx_channel,:);
                    if para_b1.dimen==3 && ~sense_philips
                        error('*** 3D SENSE ref scan is not supported! ***');
                    end
                    k_pimg_tmp=fftmr(k_pimg_tmp,1,[1:para_b1.dimen]);
                    para_b1.n_channels=length(idx_channel);
                end
                field_name=['mid',...
                            int2str(mid_pimg(imid_pimg))];
                if imid_pimg==1
                    k_pimg=struct(field_name,k_pimg_tmp);
                    para_pimg=struct(field_name,para_b1);
                else
                    k_pimg.(field_name)=k_pimg_tmp;
                    para_pimg.(field_name)=para_b1;
                end
            end
            s.mid_pimg=mid_pimg;
            s.k_pimg=k_pimg;
            s.para_pimg=para_pimg;
            if sense_philips
                % use sense ref from philips for recon
                % this is rare
                s.sense_philips=1;
            end
        end
        % save data
        if ~no_save
            save_data(fn_out,s);
            file_permission(fn_out,'+rw','ugo');
        end
        % prepare main acqusition
        if ~no_main
            if no_b0_main
                fn_main=rp(['mid',num2str(mid(imid)),...
                            '.k_nnav_uncomb.svd']);
            else
                fn_main=rp(['mid',num2str(mid(imid)),...
                            '.k_nav0_uncomb.svd']);
            end
            if strcmp(vendor,'siemens')
                disp('*** Processing the main GRE/EPI data ***');
                if ~isempty(mid_blpo)
                    kmain=recon_amri_epi(mid(imid),'k_return',1,'mid_blpo',mid_blpo,...
                                         'no_b0_crct',no_b0_main,'ste',s);
                else
                    kmain=recon_amri_epi(mid(imid),'k_return',1,...
                                         'no_b0_crct',no_b0_main,'ste',s);
                end
                kmain=kmain(:,:,:,:,find(mask_channel),:);
            elseif strcmp(vendor,'philips')
                % recon_bmir_epi doesn't perform B0 correction
                kmain=recon_bmir_epi(mid(imid),'k_return',1);
                kmain=kmain(:,:,:,find(mask_channel),:);
                % b0 correction
                if ~no_b0_main
                    df=s.df_shot;
                    if para.isgre
                        % gre
                        te=para.te_contr(:)*1e-3;
                        dp=2*pi*te.*df(:).';
                        dp=reshape(dp,[1,length(te),1,1,length(df)]);
                    else
                        % epi, echo shifting
                        % 1. base line echo time for each line
                        te=para.te_ro(:)-para.te+para.te_contr(:).';
                        te=te(:)*1e-3;
                        % 2. consider the order how echo time shift
                        n_interl=para.n_interleaves;
                        necho=length(para.te_contr);
                        nshot=length(df);
                        if isfield(para,'pe_order') && strcmp(para.pe_order,'rev_linear')
                            % for philips, the pe order can be reversed
                            te=te+...
                               [n_interl-1:-1:0]*...
                               para.echo_spacing*1e-6/n_interl*para.int_te_shift;
                        else
                            te=te+...
                               [0:n_interl-1]*...
                               para.echo_spacing*1e-6/n_interl*para.int_te_shift;
                        end
                        n_par_rep=para.n_partitions/para.sense_rate_s*para.n_reps;
                        if isfield(para,'loop_order') && strcmp(para.loop_order,'zy_order')
                            te=repmat(te,[n_par_rep,1]);
                        else
                            te=repmat(te,[1,n_par_rep]);
                        end
                        te=reshape(te,[1,para.nk_shot*necho,1,1,nshot]);
                        dp=2*pi*te.*reshape(df,[1,1,1,1,nshot]);
                    end
                    kmain=kmain./exp(1i*dp);
                end
            end
            save_data(fn_main,kmain);
            file_permission(fn_main,'+rw','ugo');
        end

    end
end
