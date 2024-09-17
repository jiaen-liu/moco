function kyz=get_pe_mdh(mid,para,acq_type)
% acq_type: 1-main, 2-blipoff
    if nargin<2
        para=extract_para(mid);
    end
    if nargin<3
        % main aquisitions
        acq_type=1;
    end
    fnmdh=get_file_filter('.',['MID*',num2str(mid),'.mdh']);
    nslice=para.n_slices;
    nrep=para.n_reps;
    nnav=0;
    if para.b_nav_en
        nnav = conditional(para.nav_type==0,1,para.nav_type);
    end
    necho_ste=0;
    if para.b_ste_en
        necho_ste=para.n_echo_steref;
    end
    necho_main=0;
    for icontr=1:para.n_contrasts
        necho_main=necho_main+...
            (para.np/para.n_interleaves/para.sense_rate_p+para.n_refs(icontr));
    end
    n_shot_noise=0;
    if para.n_noise_shots>0
        n_shot_noise=para.n_noise_shots*nslice;
    end
    n_shot_blipoff=0;
    if para.n_blipoff_reps>0
        n_shot_blipoff=para.n_blipoff_tr*nslice;
    end
    
    if acq_type==1
        n_shot_pre=n_shot_blipoff+n_shot_noise;
        nrep=para.n_reps;
        nshot_slab=para.n_partitions/para.sense_rate_s*...
            para.n_interleaves;
        if ~para.ir_bunch_slices || ~para.bunch_no_cyc
            nshot_slab=nshot_slab*para.n_delays;
        end
        nshot=nshot_slab*nslice*nrep+n_shot_blipoff+...
          n_shot_noise;
    elseif acq_type==2
        n_shot_pre=n_shot_noise;
        if para.time_stamp_date>220411
            nrep=1;
        else
            nrep=para.n_blipoff_reps;
        end
        nshot_slab=n_shot_blipoff;
        nshot=n_shot_blipoff+n_shot_noise;
    end

    
    nch=para.n_channels;
    nechot = nnav+necho_ste+necho_main;
    
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
    % get ky and kz for main acquisition
    kyz=zeros(2,necho_main,nslice,nshot_slab,nrep);
    nshot_vol=nslice*nshot_slab;

    for ir=1:nrep
        for ishot=1:nshot_slab
            for is=1:nslice
                for ie=1:necho_main
                    mdhtmp=cast2struct(dmdh(:,1,nnav+necho_ste+ie,...
                                            is+(ishot-1)*nslice+(ir-1)*nshot_vol+...
                                            n_shot_pre),...
                                       mdh);
                    kyz(1,ie,is,ishot,ir)=mdhtmp.ky;
                    kyz(2,ie,is,ishot,ir)=mdhtmp.kz;
                end
            end
        end
    end
end
