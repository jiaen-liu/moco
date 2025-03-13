function par=ste_default_conf
% par.isen=100;
    par.n_seg_mem=1;
    par.poly_fit=0;
    par.poly_fit_ker=20;
    par.cp_z_shift=0;
    par.frac_sense_sm=8;
    par.save_temp=0;
    par.simul_noise=0;
    par.npca=1;
    par.pcaflag=0;
    par.max_nch_seg=32;
    par.max_f=3;
    par.nc_max=6;
    par.ord=6;
    par.nb0=0;
    par.gb0=1;
    par.n_gb0_db0=0;
    par.ncorrection=0;
    par.uw=1;
    par.fast_uw=1;
    par.use_pri_b0=0;
    par.use_pri_b1=0;
    par.use_pri_mot=0;
    par.crct_eddy=0;
    par.mecho=0;
    par.mot_cmb=1; % 1: combine, 0: fast navigator; 2: full navigator
    par.shift_ste=1;
    % time constant for low-pass filter of motion
    % avoid filtering out respiration which is around 3 s
    par.t_mot_filt=2.5; 
    par.db0_f_filter=[0.5,0.5,0.5,0.0001];
    % clustering based on b0 maps
    par.cluster_b0=1;
    % automatically determine cluster number
    par.cluster_b0_auto=0.1;
    % pick the reference volume 
    par.auto_isen=1;
    par.sense_necho=1;
    %% less likely changed
    par.rect_cyc=1; % 0: average, 1: regression
    par.max_t=0.2;
    par.max_a=0.2;
    par.max_a=par.max_a/180*pi;% rad
    par.max_t_km=1;
    par.max_a_km=1;
    par.max_a_km=par.max_a_km/180*pi;% rad
    par.mask_thrd_ste=0.3;
    par.mask_thrd_ste_p=8e-6;
    par.nc0=1;
    par.mot_method='amri';
    par.n_iter_reset=3;
    par.n_iter=2;
    par.kos=2;
    par.en_parfor=0;
    par.ste_intp_res=4;
    par.precision='single';
    par.discard=0;
    par.nkx_g=3;
    par.nky_g=2;
    par.nkz_g=2;
    if par.use_pri_mot==1
        par.mot_cmb=0;
    end
    if par.use_pri_b0==1
        par.n_gb0_db0=1;
    end
end
