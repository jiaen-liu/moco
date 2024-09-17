function para=extract_para(mid,siemens_only,dir)
    vendor='siemens';
    if nargin<2
        siemens_only=0;
    end
    if nargin<3
        dir=pwd;
    end
    dir_cur=pwd;
    cd(dir);
    if ~isstruct(mid)
        header=siemens_asc_header(mid);
    else
        header=mid;
    end
    freq = header.sTXSPEC.asNucleusInfo{1}.lFrequency;
    frequency=freq;
    idea_v = hr_ideaversion(header.ulVersion);
    % number of contrasts
    n_contrasts = header.lContrasts;
    fa = header.adFlipAngleDegree{1};
    tr = header.alTR{1}/1000.0;
    if isfield(header, 'alTI') 
        ti=zeros(length(header.alTI),1);
        for i=1:length(header.alTI)
            ti(i) = header.alTI{i}/1000.0;
        end
    else
        ti = [];
    end;
    te = header.alTE{1}/1000.0;
    te_contr=zeros(n_contrasts,1);
    for i=1:n_contrasts
        te_contr(i)=header.alTE{i}/1000;
    end
    nte_contr=length(te_contr);
    b_fatsat = false;
    if header.sPrepPulses.ucFatSat==1 
        b_fatsat = true;
    end
    n_channels = 0;
    if contains(idea_v, 'VB', 'IgnoreCase',true)
        idx_channels=zeros(length(header.asCoilSelectMeas{1}.asList),1);
        for i=1:length(header.asCoilSelectMeas{1}.asList)
            if header.asCoilSelectMeas{1}.asList{i}.lElementSelected
                n_channels=n_channels+1;
                idx_channels(i)=1;
            else
                idx_channels(i)=0;
            end
        end
    elseif contains(idea_v, 'VD', 'IgnoreCase',true) || ...
            contains(idea_v, 'VE', 'IgnoreCase',true)
        idx_channels=zeros(length(header.sCoilSelectMeas.aRxCoilSelectData{1}.asList),1);
        name_channels=string;
        for i=1:length(header.sCoilSelectMeas.aRxCoilSelectData{1}.asList)
            if header.sCoilSelectMeas.aRxCoilSelectData{1}.asList{i}.lElementSelected
                n_channels=n_channels+1;
                idx_channels(i)=header.sCoilSelectMeas.aRxCoilSelectData{1}.asList{i}.lADCChannelConnected;
                name_channels(i)=header.sCoilSelectMeas.aRxCoilSelectData{1}.asList{i}.sCoilElementID.tElement;
                % idx_channels(i)=header.sCoilSelectMeas.aRxCoilSelectData{1}.asList{i}.sCoilElementID.ulUniqueKey;
            else
                idx_channels(i)=0;
            end
        end
    end
    multislicemode = header.sKSpace.ucMultiSliceMode;
    slice_order = header.sSliceArray.ucMode;
    % number of repetitions
    n_reps = 1;
    if isfield(header, 'lRepetitions')
        n_reps = header.lRepetitions+1;
    end
    % number of averages
    n_aves = 1;
    if isfield(header, 'lAverages') 
        n_aves = header.lAverages;
    end
    % number of reptions is just multipled
    % by number of averages
    % dwell time (time between samples)
    t_dwell = header.sRXSPEC.alDwellTime{1};
    nr_os = 2.0*header.sKSpace.lBaseResolution;
    nr = header.sKSpace.lBaseResolution;
    % phase lines
    np = header.sKSpace.lPhaseEncodingLines;
    % nuber of slices (2d) is the same as slabs (3d)
    n_slices = header.sSliceArray.lSize;
    sDistFactor = 0;
    if isfield(header.sGroupArray.asGroup{1}, 'dDistFact')
        sDistFactor = header.sGroupArray.asGroup{1}.dDistFact;
    end
    if n_slices == 1 
        sDistFactor = 0;
    end
    sl_oversamp = 0;
    if isfield(header.sKSpace, 'dSliceOversamplingForDialog')
        sl_oversamp = header.sKSpace.dSliceOversamplingForDialog;
    end
    % number of partitions
    n_partitions = 1;
    n_partitions_nos = 1;
    dimen = 2;
    if header.sKSpace.ucDimension == 4
        n_partitions = round(header.sKSpace.lImagesPerSlab*(1+sl_oversamp));
        n_partitions_nos = header.sKSpace.lImagesPerSlab;
        dimen = 3;
    end
    % fov in mm
    fovr = header.sSliceArray.asSlice{1}.dReadoutFOV;
    fovp = header.sSliceArray.asSlice{1}.dPhaseFOV;
    sthickness = header.sSliceArray.asSlice{1}.dThickness;
    resr = double(fovr)/nr;
    resp = double(fovp)/np;
    ress = double(sthickness)/n_partitions_nos*(1+sDistFactor);
    % rotation
    if isfield(header.sSliceArray.asSlice{1},'dInPlaneRot')
        prot=header.sSliceArray.asSlice{1}.dInPlaneRot;
    else
        prot=0;
    end
    % fov shift
    x_shift = zeros(n_slices,1);
    y_shift = zeros(n_slices,1);
    z_shift = zeros(n_slices,1);
    r_shift = zeros(n_slices,1);
    p_shift = zeros(n_slices,1);
    s_shift = zeros(n_slices,1);
    % tra, sag, cor
    snormal = zeros(3, n_slices);
    for i = 1:n_slices
        if isfield(header.sSliceArray.asSlice{i}, 'sPosition')
            if isfield(header.sSliceArray.asSlice{i}.sPosition, 'dTra') 
                z_shift(i) = ...
                    header.sSliceArray.asSlice{i}.sPosition.dTra;
            end
            if isfield(header.sSliceArray.asSlice{i}.sPosition, 'dSag') 
                x_shift(i) = ...
                    header.sSliceArray.asSlice{i}.sPosition.dSag;
            end
            if isfield(header.sSliceArray.asSlice{i}.sPosition, 'dCor') 
                y_shift(i) = ...
                    header.sSliceArray.asSlice{i}.sPosition.dCor;
            end
        end

        if isfield(header.sSliceArray.asSlice{i}.sNormal, 'dSag') 
            snormal(1, i) = ...
                header.sSliceArray.asSlice{i}.sNormal.dSag;
        end
        if isfield(header.sSliceArray.asSlice{i}.sNormal, 'dCor') 
            snormal(2, i) = ...
                header.sSliceArray.asSlice{i}.sNormal.dCor;
        end
        if isfield(header.sSliceArray.asSlice{i}.sNormal, 'dTra') 
            snormal(3, i) = ...
                header.sSliceArray.asSlice{i}.sNormal.dTra;
        end
    end
    % to do r, p and s shift should not be equal to x, y and z
    r_shift = x_shift;
    p_shift = y_shift;
    s_shift = z_shift;
    para=var2struct('freq','frequency','idea_v','n_contrasts',...
                    'fa','tr','ti','te','te_contr','nte_contr',...
                    'b_fatsat','n_channels','multislicemode',...
                    'n_reps','n_aves','t_dwell',...
                    'nr_os','nr','np','n_slices',...
                    'sDistFactor','sl_oversamp','n_partitions_nos',...
                    'n_partitions','dimen',...
                    'fovr','fovp','sthickness','snormal',...
                    'resr','resp','ress','prot',...
                    'x_shift','y_shift','z_shift',...
                    'r_shift','p_shift','s_shift','vendor','slice_order');
    if ~siemens_only
        wip=amri_epi_wipmem(header);
        wip=struct2double(wip);
        t_dwell_ste = 1.0/wip.steref_bandwidth*1e9/2;
        % inversion pulse train length
        ir_pulse_train_length = wip.ir_pulse_train_length;
        % ms
        ir_duration = wip.ir_duration_us/1000.0; 
        % mt pulse flip angle
        mtc_fa = wip.mtc_angle_deg;
        % rf time bandwidth product
        rf_tbw = wip.rf_bandw_time_prod;
        % fat saturation pulse flip angle
        fatsat_fa = wip.fatsat_angle_deg;
        % 
        time_stamp = wip.compilation_timestamp;
        time_stamp_date=floor(time_stamp/10000);
        % the unit in wip was 10 us
        vte_arr = [0;wip.var_te_array]/100; 
        % navigator
        b_nav_en = wip.navigator;
        n_navs = zeros(12,1,'int16');
        b_nav_pol = zeros(12,1,'int16');
        nav_type = wip.nav_type;
        nav_bw = wip.nav_bandwidth;
        nav_dim = [wip.nav_dim_r;wip.nav_dim_p;wip.nav_dim_s];
        % epi polarity
        b_epi_pol = zeros(12,1,'int16');
        % epi positive only readout
        b_epi_positive = wip.positive_only;
        fixed_n_interleaves = wip.fixed_n_interleaves;
        for i = 1:n_contrasts
            if b_nav_en 
                n_navs(i) = bitand(bitshift(wip.navigator_bitmask, -i+1),int32(1));
            end
            b_nav_pol(i) = bitand(bitshift(wip.nav_polarity_bitmask, -i+1),int32(1));
            b_epi_pol(i) = bitand(bitshift(wip.epi_polarity_bitmask, -i+1),int32(1));
        end
        n_echo_nav=0;
        if b_nav_en
            n_echo_nav = conditional(nav_type==0,1,nav_type);
        end
        % number of epi interleaves
        n_interleaves = wip.n_interleaves;
        % sense
        sense_rate_p = wip.sense_rate_p;
        sense_rate_s = wip.sense_rate_s;
        % bydefault, this is zero unless
        % defined
        dkz_caipi = 0;
        if wip.caipi==1
            dkz_caipi=floor(sense_rate_s/2);
        end
        n_sense_reps = wip.n_sense_ref;
        % number of variable TEs
        n_var_tes = wip.n_var_te;
        n_sense_tr = n_sense_reps*n_partitions*...
            n_interleaves*sense_rate_p*n_var_tes;
        % indicate the state of cyclic variable te
        b_cyc_var_te = wip.var_te_cyclic;
        % indicate the state of fast/slow variable te loop
        b_fast_var_te = wip.var_te_fastloop;
        n_fvarte_cyc = 1;
        n_varte_cyc = 1;
        n_fvarte_ncyc = 0;
        n_varte_ncyc = 0;
        % determine the respective number of loops for variable TEs
        if n_var_tes > 1 && b_cyc_var_te
            if b_fast_var_te 
                n_fvarte_cyc = n_var_tes;
            else
                n_varte_cyc = n_var_tes;
            end

        elseif n_var_tes > 1
            if b_fast_var_te
                n_fvarte_ncyc = n_var_tes;
                n_varte_ncyc = 1;
            else
                n_fvarte_ncyc = 1;
                n_varte_ncyc = n_var_tes;
            end
        end
        % non-cyclic te loops
        n_nvar_tes = n_var_tes;
        n_nvarte_vol = 0;
        if ~wip.var_te_cyclic && wip.var_te
            n_nvarte_vol = 1;
        end
        n_nvarte_tr = ...
            n_nvarte_vol*n_partitions/sense_rate_s*...
            n_interleaves*n_nvar_tes;
        % ramp up or down time
        ramp_dur = wip.ramp_duration_us;
        % ramp sampling fraction
        if isfield(wip,'rampsampling')
            % new sequence modified on May 2021
            max_necho=12;
            if isempty(header.alTE{127-5*max_necho+1})
                ramp_samp_frac=0;
            else
                ramp_samp_frac=double(header.alTE{127-5*max_necho+1})/1e6;
            end
        else
            ramp_samp_frac = wip.rampsamp_frac;
        end
        if n_interleaves*sense_rate_p == np 
            isgre = 1;
        else 
            isgre = 0;
        end
        % number of inversion/recovery delays
        n_delays = wip.ir_n_ti;
        % delay times
        ti = ti(1:n_delays);
        % number of blip off repetitions
        n_blipoff_reps = wip.n_blipoff_vol;
        if time_stamp_date>220411
            n_blipoff_tr = min(n_blipoff_reps,...
                               n_partitions/sense_rate_s);
        else
            n_blipoff_tr = n_blipoff_reps*...
                n_partitions/sense_rate_s;
        end
        % number of reference lines per contrast
        n_refs = wip.n_ref_echoes;

        % number of k-lines per shot
        nk_shot = np/n_interleaves/sense_rate_p;
        nk_shot_ref = nk_shot+n_refs(1);
        % fixed interleaevs
        b_fix_n_interleaves = wip.fixed_n_interleaves;
        % time between echos
        echo_spacing = wip.echo_spacing;
        int_te_shift = wip.int_te_shift;
        
        if isgre
            nte_contr=0;
            te_contr=[];
            for i=1:n_contrasts
                nte_contr=nte_contr+n_refs(i)+1;
                te_contr=[te_contr;...
                          header.alTE{i}/1000+...
                          [0:n_refs(i)].'*echo_spacing*1e-3];
            end
        end        

        if isfield(wip,'ste3d_mode')
            ste3d_mode = wip.ste3d_mode;
        else 
            ste3d_mode = 0;
        end
        ste_rsamp = wip.ste_rsamp/10000;
        steref_dim_r = wip.steref_dim_r;
        steref_dim_p = wip.steref_dim_p;
        steref_dim_s = wip.steref_dim_s;
        steref_res_r = fovr/steref_dim_r;
        steref_res_p = fovp/steref_dim_p;
        steref_res_s = sthickness*(1.0+sl_oversamp)/steref_dim_s;
        n_interleaves_steref = wip.n_interleaves_steref;
        n_echo_steref = 0;
        b_ste_en = 0;
        if wip.short_te_ref
            n_echo_steref = steref_dim_p/n_interleaves_steref;
            b_ste_en = 1;
        end
        dte_ste=0;
        te_ste=0;
        if b_ste_en
            dte_ste=wip.ste_rtime*1e-3+1/wip.steref_bandwidth*steref_dim_r*1000;
            te_ste=double(wip.echo_time_steref_us)*1e-3;
            te_ste=([0:n_echo_steref-1]-floor(n_echo_steref/2))*dte_ste+te_ste;
        end
        idx_nref = (1:floor(nk_shot/2)+1).';
        if floor(nk_shot/2)+n_refs(1)+1 <= nk_shot_ref-1
            idx_nref = [idx_nref;[floor(nk_shot/2)+n_refs(1)+2:nk_shot_ref].'];
        end
        % polarity of each k-space line within one shot
        % 1 means negative, 0 positive
        ro_pol_ref = zeros(nk_shot_ref, n_contrasts);
        para_tmp = struct('b_epi_positive',b_epi_positive,...
                          'b_epi_pol',b_epi_pol);
        for i=1:nk_shot_ref
            for j=1:n_contrasts
                ro_pol_ref(i,j)=line_pol(i,j,para_tmp);
            end
        end
        ro_pol = ro_pol_ref(idx_nref, :);
        % in millisecond
        te_ro_ref = ((0:nk_shot_ref-1).'-floor(nk_shot/2))*...
            echo_spacing*1e-3+te;
        te_ro = te_ro_ref(idx_nref);
        % count TRs
        ir_bunch_slices=wip.ir_bunch_slices;
        bunch_no_cyc=wip.bunch_no_cyc;
        if wip.ir_bunch_slices && wip.bunch_no_cyc
            n_main_shots = ...
                n_reps*n_partitions/sense_rate_s*...
                n_interleaves*n_slices*(n_fvarte_cyc*n_varte_cyc);
        else
            n_main_shots = ...
                n_reps*n_delays*n_partitions/sense_rate_s*...
                n_interleaves*n_slices*(n_fvarte_cyc*n_varte_cyc);
        end
        n_main_tr = n_main_shots/n_slices;
        n_noise_tr = wip.n_noiserep;
        n_noise_shots = n_noise_tr;
        n_dummy_tr = wip.dummy_shots;
        n_dummy_shots = n_dummy_tr*n_slices;
        n_tr = n_noise_tr+n_dummy_tr+...
               n_blipoff_tr+n_sense_tr+...
               n_nvarte_tr+n_main_tr;
        version=wip.seq_version;
        parawip=var2struct('t_dwell_ste','ir_pulse_train_length','ir_duration',...
                           'mtc_fa','rf_tbw','fatsat_fa',...
                           'time_stamp','time_stamp_date','vte_arr',...
                           'b_nav_en','n_navs','b_nav_pol',...
                           'nav_type','nav_bw','nav_dim',...
                           'b_epi_pol','b_epi_positive','fixed_n_interleaves',...
                           'n_navs','b_nav_pol','b_epi_pol',...
                           'n_interleaves','sense_rate_p','sense_rate_s',...
                           'dkz_caipi','n_sense_reps',...
                           'n_var_tes','n_sense_tr',...
                           'b_cyc_var_te','b_fast_var_te',...
                           'n_fvarte_cyc','n_varte_cyc','n_fvarte_ncyc','n_varte_ncyc',...
                           'n_nvar_tes','n_nvarte_vol','n_nvarte_tr',...
                           'ramp_dur','ramp_samp_frac',...
                           'isgre','n_delays','ti',...
                           'n_blipoff_reps','n_blipoff_tr','n_refs',...
                           'nk_shot','nk_shot_ref',...
                           'b_fix_n_interleaves',...
                           'echo_spacing','int_te_shift',...
                           'ste3d_mode','ste_rsamp',...
                           'steref_dim_r','steref_dim_p','steref_dim_s',...
                           'steref_res_r','steref_res_p','steref_res_s',...
                           'n_interleaves_steref','n_echo_steref','b_ste_en',...
                           'dte_ste','te_ste',...
                           'ro_pol_ref','ro_pol',...
                           'te_ro_ref','te_ro',...
                           'n_main_tr','n_noise_tr','n_noise_shots',...
                           'n_dummy_tr','n_dummy_shots',...
                           'n_tr','nte_contr','te_contr','n_main_shots',...
                           'ir_bunch_slices', 'bunch_no_cyc','version');
        para=pass_var_struct(para,parawip);
    end
    para=struct2double(para);
    cd(dir_cur);
end
