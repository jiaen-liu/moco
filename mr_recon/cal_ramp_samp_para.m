function [rsamp,nk,flat,ramp]=cal_ramp_samp_para(para,ste)
% $$$ ; Returns (in calling vars) actual rampsample fraction, number acq. of points
% $$$ ;      flat time and ramp time (in units of dwell time)
% $$$ ;  based on siemens header and calculations in AMRI_readout
% $$$ ; Input: siemens header from siemens_amri_epi (needs both ascii and wip parts)
% $$$ ;        echo number (for selection of bandwidth from header)
% $$$ ;
    if nargin<2
        ste=false;
    end
    if ste
        dt=para.t_dwell_ste/1000.0;
        ramp=para.ramp_dur/dt;
        rsamp=para.ste_rsamp;
        nr=para.steref_dim_r*2;
    else
        dt = para.t_dwell/1000.0; % recv samp time is in ns
        ramp = para.ramp_dur/dt;    % ramp time in sample units (recv samp time is in ns)
        rsamp = para.ramp_samp_frac; % requested ramp fraction
        nr = para.nr_os; % desired res, with oversampling
    end
    flat10= (nr*dt+ ramp*(rsamp^2- 2* rsamp)*dt)/10.0;

    if (flat10 > floor(flat10)) 
        flat= 10*floor(flat10)+10;
    else 
        flat= 10*floor(flat10);
    end
    flat= flat/dt;
    if flat >= nr
        nk=nr;
        rsamp=0.0;
        return;
    end
    rsamp= 1- (1- (nr- flat)/ramp)^0.5; % actual ramp fraction
    nk= flat+ 2.0*ramp*rsamp;
    if (nk > floor(nk)) 
        nk= floor(nk)+ 1;
    else 
        nk= floor(nk);
    end
    if mod(nk,4) > 0
        nk= nk+4- mod(nk, 4); % rounded to integer * 4 for oversampling 2
    end
    rsamp= (nk-flat)*0.5/ramp;
end
