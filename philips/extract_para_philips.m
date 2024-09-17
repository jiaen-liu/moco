function para=extract_para_philips(mid)
if ischar(mid)
    fn_raw=mid;
elseif isnumeric(mid)
    fn_raw=mid2filename_philips(mid);
end
vendor='philips';
% use reconframe
r=MRecon(fn_raw);
% p=ParameterReader(fn_raw);
% $$$ r.SortData;
nkr=r.Parameter.Encoding.KxRange(1,2)-r.Parameter.Encoding.KxRange(1,1)+1;
nkp=r.Parameter.Encoding.KyRange(1,2)-r.Parameter.Encoding.KyRange(1,1)+1;
if isempty(r.Parameter.Encoding.KzRange)
    nks=1;
else
    nks=r.Parameter.Encoding.KzRange(1,2)-r.Parameter.Encoding.KzRange(1,1)+1;
end
tfe_factor=1;
if strcmp(r.Parameter.GetValue('VAL01_ACQ_tfe_excitations'),'yes')
    tfe_factor=r.Parameter.GetValue('RC_turbo_factor');
end
n_delays=r.Parameter.GetValue('VAL01_CARD_user_def_phases');

try
    freq=r.Parameter.GetValue('PR_FF_found_F0');
catch me
    if strcmp(me.message,'The parameter does not exist')
        warning(['*** PR_FF_found_F0 does not exist ***']);
        freq=[];
    end
end
frequency=freq;
try
    b_epi_positive=strcmp(r.Parameter.GetValue('VAL01_ACQ_multiecho_flyback'),'yes');
catch me
    if strcmp(me.message,'The parameter does not exist')
        warning(['*** VAL01_ACQ_multiecho_flyback does not exist ***']);
        b_epi_positive=1;
    end
end
b_epi_pol=zeros(12,1);
b_epi_pol(2:2:end)=1;
fat_shift_dir=r.Parameter.GetValue('LCA::ima->fat_shift_dir.cmp(1)');
isgre=strcmp(r.Parameter.GetValue('VAL01_ACQ_epi_readout'),'no');
sampled_fovs=r.Parameter.GetValue('LCA::ima->sampled_fovs.cmp(1)');
sense_factors=r.Parameter.GetValue('LCA::ima->sense_factors.cmp(1)');
oversample_factors=r.Parameter.GetValue('LCA::ima->oversample_factors.cmp(1)');
dimen=r.Parameter.GetValue('ENC::ima->nr_encodings.cmp(1)');
if dimen==3
    sDistFactor=0;
else
    sDistFactor=r.Parameter.GetValue('STACK::ima->slice_gap.cmp(1)');
end
n_slices=r.Parameter.GetValue('LCA::ima->comp_elements.cmp(1)');
fa=r.Parameter.GetValue('VAL01_ACQ_flip_angle');
tr=r.Parameter.GetValue('IF_act_FFE_rep_time');
TR=r.Parameter.GetValue('VAL01_CARD_effective_TR');
te=r.Parameter.GetValue('VAL05_FFE_act_first_echo_time');
grad_ro_str=r.Parameter.GetValue('GR::m_0_->str');
nte_contr=r.Parameter.GetValue('VAL01_ACQ_echoes');
channel_id=double(r.Parameter.Labels.CoilNrs(:,1));
n_channels=length(channel_id);
% r.Parameter.GetValue('AQ::base->nr_channels_enabled.cmp(1)');
% scanner coordinate: x (P->A), y (R->L), z (H->F)
omatrix_xyz=r.Parameter.GetValue('O_MATRIX::locations_0_->matrix');
omatrix_xyz=reshape(omatrix_xyz,[3,3]);
% patient coordinate: R->L, A->P, F->H
omatrix_pat=inv([0,-1,0;1,0,0;0,0,-1])*omatrix_xyz;
snormal=repmat(omatrix_pat(:),[1,n_slices]);
% prot=zeros(1,n_slices);
prot=0;
% offcenter, defined in patient coordinate
% because x direction is up-down in the Philps definition
% using patient coordinate is consistent with Siemens
x_shift=zeros(n_slices,1);
y_shift=zeros(n_slices,1);
z_shift=zeros(n_slices,1);
m_shift=zeros(n_slices,1);
p_shift=zeros(n_slices,1);
s_shift=zeros(n_slices,1);
m_freq_factor=zeros(n_slices,1);
if n_slices==1
    v=r.Parameter.GetValue('LCA::ima->pat_offcentres');
    y_shift=v(1);
    x_shift=v(2);
    z_shift=v(3);    
    v=r.Parameter.GetValue('LCA::ima->mps_offcentres');
    m_shift=v(1);
    p_shift=v(2);
    s_shift=v(3);
    m_freq_factor=r.Parameter.GetValue('LCA::ima->lca_m_freq_factor');
else
    for i=1:n_slices
        v=r.Parameter.GetValue(['LCA::ima->pat_offcentres.cmp(',num2str(i),')']);
        y_shift(i)=v(1);
        x_shift(i)=v(2);
        z_shift(i)=v(3);
        v=r.Parameter.GetValue(['LCA::ima->mps_offcentres.cmp(',num2str(i),')']);
        m_shift(i)=v(1);
        p_shift(i)=v(2);
        s_shift(i)=v(3);
        m_freq_factor(i)=r.Parameter.GetValue(['LCA::ima->lca_m_freq_factor.cmp(',num2str(i),')']);
    end
end
clear v i;
% ns
t_dwell=r.Parameter.GetValue('AQ::base->interval.cmp(1)')*1e6;
n_dyn=r.Parameter.GetValue('VAL01_DYN_nr_scans');
n_meas=r.Parameter.GetValue('VAL01_ACQ_measurements');
n_reps=n_dyn*n_meas;
n_aves=n_meas;
n_refs=zeros(12,1);
if isgre
    n_refs(1)=nte_contr-1;
    n_contrasts=1;
else
    n_contrasts=nte_contr;
end
epi_factor=1;
if strcmp(r.Parameter.GetValue('VAL01_ACQ_fast_imaging_mode'),'EPI') || strcmp(r.Parameter.GetValue('VAL01_ACQ_fast_imaging_mode'),'TFEPI')
   if ~strcmp(r.Parameter.GetValue('VAL01_ACQ_shot_mode'),'single-shot')
       epi_factor=r.Parameter.GetValue('VAL01_ACQ_epi_factor');
   else
       epi_factor=nkp;
   end
end
nk_shot=epi_factor;
nk_shot_ref=nk_shot+n_refs(1);
int_te_shift=strcmp(r.Parameter.GetValue('CSC_epi_echo_shift'),'yes');
pe_order=r.Parameter.GetValue('CSC_ffe_profile_order_sel');
loop_order=r.Parameter.GetValue('VAL01_FFE_loop_order');
%%
dte=r.Parameter.GetValue('VAL05_FFE_echo_time_step');
if isgre
    echo_spacing=dte*1000;
else
    echo_spacing=r.Parameter.GetValue('GR::m_0_->dur')*1000;
end
te_contr=[0:nte_contr-1]*dte+te;
sl_oversamp=max(0,oversample_factors(3)-1);
sense_factors_int=round(max(1,sense_factors./oversample_factors));
sense_rate_p=sense_factors_int(2);
% sense_rate_s=sense_factors_int(3);
if mod(sense_factors(3),1)~=0
    % error('*** This scan contains extra oversampling in the slice directin! ***');
    sense_rate_s=round(sense_factors(3));
    sense_factors_int(3)=sense_rate_s;
else
    sense_rate_s=sense_factors(3);
    sense_factors_int(3)=sense_factors(3);
end
dkz_caipi=0;
% mm unit
resr=sampled_fovs(1)/nkr;
resp=sampled_fovs(2)/nkp;
ress=sampled_fovs(3)/nks*(1+sDistFactor);
fov=sampled_fovs.*sense_factors_int;
fov(1)=fov(1)/oversample_factors(1);
fovr=fov(1);
fovp=fov(2);
fovs=fov(3);
sthickness=fov(3)/(1+sl_oversamp);
nr=round(fovr/resr);
nr_os=nr*oversample_factors(1);
nsamp_ro=r.Parameter.GetValue('AQ::base->samples.cmp(1)');
np=round(fovp/resp);
n_partitions=nks*sense_rate_s;
n_partitions_nos=round(n_partitions/(1+sl_oversamp));
if isgre
    n_interleaves=np/sense_rate_p;
else
    n_interleaves=np/epi_factor/sense_rate_p;
end
nus_enc=r.Parameter.GetValue('RC_nus_enc_nrs');
nus_enc=nus_enc(1:min(nsamp_ro,length(nus_enc)));
ramp_dur=r.Parameter.GetValue('GR::m_0_->slope')*1000;
ramp_samp_frac=r.Parameter.GetValue('VAL02_ACQ_epi_aq_slope_fraction');
if ramp_samp_frac>1e-3 && all(nus_enc==0)
    ramp_samp_frac=0;
end
% rx_channel_id=r.Parameter.Labels.CoilNrs(:,1);
%%
b_nav_en=0;
n_navs=zeros(12,1);
b_nav_pol=zeros(12,1);
nav_type=0;
nav_bw=0;
nav_dim=[0,0,0];
n_sense_reps=0;
n_sense_tr=0;
n_varte_vol=0;
n_varte_tr=0;
n_blipoff_reps=0;
n_blipoff_tr=0;
n_noise_tr=0;
n_noise_shots=0;
n_dummy_tr=0;
n_dummy_shots=0;
multislicemode=1;
b_fatsatsat=0;
%%
% nav3d
b_ste_en=0;
try
    b_ste_en=~(strcmp(r.Parameter.GetValue('VAL01_FFE_nav3d'),'none') || ...
               strcmp(r.Parameter.GetValue('VAL01_FFE_nav3d'),'MGU_NAV3D_NONE'));
catch
    b_ste_en=0;
end
if b_ste_en
    steref_dim_r=r.Parameter.GetValue('VAL01_FFE_nr_nav3d');
    steref_dim_p=r.Parameter.GetValue('VAL01_FFE_np_nav3d');
    steref_dim_s=r.Parameter.GetValue('VAL01_FFE_ns_nav3d');
    steref_s1=r.Parameter.GetValue('VAL01_FFE_sy_nav3d');
    steref_s2=r.Parameter.GetValue('VAL01_FFE_sz_nav3d');
    ste_rsamp=0;
    ste3d_mode=101;
    steref_res_r=fovr/steref_dim_r;
    steref_res_p=fovp/steref_dim_p;
    steref_res_s=sthickness*(1+sl_oversamp)/steref_dim_s;
    n_interleaves_steref=r.Parameter.GetValue('VAL01_FFE_sy_nav3d');
    n_echo_steref=steref_dim_p/n_interleaves_steref;
    t_dwell_ste=r.Parameter.GetValue('AQ::nav3d->interval.cmp(1)')*1e6;
    te_ste=r.Parameter.GetValue('VAL05_FFE_act_nav3d_time');
    dte_ste=r.Parameter.GetValue('GR::m_nav3d->dur.cmp(1)');
    te_ste=([0:n_echo_steref-1]-floor(n_echo_steref/2))*dte_ste+te_ste;
end
%%
te_ro=((0:nk_shot-1).'-floor(nk_shot/2))*...
          echo_spacing*1e-3+te;
% te_ro and te_ro_ref should only used for EPI data
te_ro_ref=te_ro;
ro_pol_ref=zeros(nk_shot_ref,n_contrasts);
if ~b_epi_positive
    for icontr=1:n_contrasts
        ro_pol_ref((2-mod(icontr-1,2)):2:end,icontr)=1;
    end
end
ro_pol=ro_pol_ref;
n_main_tr=r.Parameter.GetValue('RC_total_nr_profiles')/...
          nte_contr/epi_factor/n_channels/n_slices;
n_tr=n_main_tr;
% $$$ channel_id=[];
% $$$ for i=1:length(r.Parameter.Labels.CoilInfo)
% $$$     channel_id=[channel_id;...
% $$$                 col(r.Parameter.Labels.CoilInfo(i).channel_uids(find(...
% $$$                     r.Parameter.Labels.CoilInfo(i).active_channels)))];
% $$$ end

% channel_id=double(channel_id);

ecc_filter_counts=r.Parameter.GetValue('HW_ecc_filter_counts');
rf_b0c_type=r.Parameter.GetValue('RF::ex->b0c_type.cmp(1)');
nav3d_freq_ref=r.Parameter.GetValue('AQ::nav3d->freq_ref.cmp(1)');
rf_ex_shape=r.Parameter.GetValue('RF::ex->shape.cmp(1)');
rf_ex_dur=r.Parameter.GetValue('RF::ex->dur.cmp(1)');
k_ellip=strcmp(r.Parameter.GetValue('VAL00_DEF_elliptical_k_space_shutter'),'yes');
rf_spoil_phase=r.Parameter.GetValue('MPF_spoil_phase_angle');

nav1d_enable=0;
te_nav1d=[];
if strcmp(r.Parameter.GetValue('VAL01_FFE_phase_nav'),'pre-acquisition')
    nav1d_enable=1;
    te_nav1d(1)=r.Parameter.GetValue('AQ::phnav_0_->time');
elseif strcmp(r.Parameter.GetValue('VAL01_FFE_phase_nav'),'post-acquisition')
    nav1d_enable=1;
    te_nav1d(1)=...
        r.Parameter.GetValue('SQ::base->dur')-...
        r.Parameter.GetValue('SQ::base->ref')+...
        r.Parameter.GetValue('SQ::xbase->dur')+...
        (nte_contr-1)*(r.Parameter.GetValue('SQ::xME->dur')+...
                       r.Parameter.GetValue('SQ::ME->dur'))+...
        r.Parameter.GetValue('SQ::fin->ref')+...
        r.Parameter.GetValue('AQ::phnav_1_->time');
elseif strcmp(r.Parameter.GetValue('VAL01_FFE_phase_nav'),'both')
    nav1d_enable=2;
    te_nav1d(1)=r.Parameter.GetValue('AQ::phnav_0_->time');
    te_nav1d(2)=...
        r.Parameter.GetValue('SQ::base->dur')-...
        r.Parameter.GetValue('SQ::base->ref')+...
        r.Parameter.GetValue('SQ::xbase->dur')+...
        (nte_contr-1)*(r.Parameter.GetValue('SQ::xME->dur')+...
                       r.Parameter.GetValue('SQ::ME->dur'))+...
        r.Parameter.GetValue('SQ::fin->ref')+...
        r.Parameter.GetValue('AQ::phnav_1_->time');
end
if nte_contr==1 && nav1d_enable==2
    te_nav1d(2)=r.Parameter.GetValue('AQ::phnav_1_->time');
end
% With 34 channels, the last two channels are tossed out for EPI phase correction specific to the 7T
if abs(frequency/42.58e6-7)<0.5
    if n_channels>32
        if n_channels~=34
            error('*** There should be 32 or 34 channels! ***');
        end
        % remove the body coil or TR coil
        channel_id=channel_id(channel_id>1);
        n_channels=length(channel_id);
    end
end
% RF drive scale
rf_ds=r.Parameter.GetValue('PR_rf_drive_scales');
% reconframe object is cleared.
clear r;
para=workspace2struct(who);


