function [noise,mag_noise]=noise_map_moco_recon(b1,cov_mat,para,para_b1,downsample,tukey_ratio,legacy)
% load b1 as load(['result/mid',num2str(mid),'.b1.mat']);
% the image result is in ['result/im',num2str(mid),'.result_MoCo_B0Co.mat']);
% cov_mat=covSiem(mid);
    if nargin<6
        tukey_ratio=[];
    end
    if nargin<7
        legacy=0;
    end
    para.resr=para.resr*downsample(1);
    para.resp=para.resp*downsample(2);
    para.ress=para.ress*downsample(3);
    
    para.nr=para.nr/downsample(1);
    para.np=para.np/downsample(2);
    para.n_partitions=para.n_partitions/downsample(3);
    para.n_partitions_nos=para.n_partitions_nos/downsample(3);
    if para.n_slices>1
        error('Multi-slice acquisition is not supported!');
    end
    s1=para.sense_rate_p;
    s2=para.sense_rate_s;
    dz=para.dkz_caipi;
    
    coord=get_coordinate(para,1,0,0);
    coord_b1=get_coordinate(para_b1,1,0,0);
    
    [nx,ny,nz,nch]=size(b1);
    b1=reshape(b1,[nx*ny*nz,nch]);
    % this is already done in ste_final_prep_nintp
    % Not needed here and commented out
    % b1=b1*conj(chol(inv(cov_mat),'lower'));
    b1=b1/prctile(abs(b1(:)),95);
    b1=reshape(b1,[nx,ny,nz,nch]);

    b1hres=interp_sense(coord_b1,coord,b1);
    [nx,ny,nz,nch]=size(b1hres);
    nx=nx*downsample(1);
    ny=ny*downsample(2);
    nz=nz*downsample(3);
    
    % sensit_ds=fft_down_sample(b1hres,[1,2,3],downsample,1,0);
    sensit_ds=b1hres*prod(downsample);
    % [~,noise_ds]=sense_mat(sensit_ds,inv(cov_mat),s1,s2,dz,0,0);
    % instead of using cov_mat, identity matrix is used becaue 
    % b1 has been normalized earlier
    [~,noise_ds]=sense_mat(sensit_ds,eye(nch,nch),s1,s2,dz,0,0);
    noise_ds=noise_ds.^0.5;
    % noise=im_intp_res(noise_ds,1./downsample)*prod(downsample)^0.5;
    
    noise=im_intp_res(noise_ds,1./downsample)*prod(downsample)^0.5;
    % because recon_epi_ste_beta use cov_mat with oversampling
    % a factor of 2^0.5 is needed
    % 10-28-2024, JL: it's actually not needed, check validate_moco_recon_noise.m
    % noise=noise/2^0.5;
    % 12-17-2024: the reason is due to the tukey fraction below which is about a factor of 1/2^0.5
    % in recon_epi_ste_beta, the simulation doesn't consider tukey filtering
    % as a result, the simulation and noise calculation in 
    % validate_moco_recon_noise.m gives similar level of noise
    % which is not correct
    % we should still use the following:
    if ~legacy
        noise=noise/2^0.5;
    end
    if ~isempty(tukey_ratio)
        noise=noise*tukey_fraction([nx,ny,nz],tukey_ratio).^0.5;
    end
    var_mag_noise=abs(noise).^2/2;
    mag_noise=var_mag_noise.^0.5;
    
end
