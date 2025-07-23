function convertBOLD(mid,TR,whichReconPar,isNegPEDir)
% Example - convertBOLD(92,2.34,'moco',0);

    if ischar(mid)
        mid=eval(mid);
    end

    switch whichReconPar
    case 'moco'
        reconName = 'MoCo_B0Co_RefB0_Sag'; %motion and B0 correction
    case 'nomoco'
        reconName = 'gB0Co_RefB0_Sag'; %B0 correction 
    case 'noco'
        reconName = 'noB0MoCo_Sag'; %no correction
    otherwise
        error('WRONG recon configuration input!');
    end
    
    fname_img_ave=get_file_filter('./result/',['im',num2str(mid(1)),'_', ...
        reconName,'_mag_echo_ave.nii.gz']); % new sort_siemens output
    fname_img_all_echoes=get_file_filter('./result/',['im',num2str(mid(1)), ...
        '.result_',reconName,'.mat']); % new sort_siemens output
    img_niiInfo = niftiinfo(fname_img_ave);
    img_all_echoes = load(fname_img_all_echoes);
    img_niiInfo_orig = img_all_echoes.par;
    
    
    img_all_echoes.im_recon = abs(squeeze(img_all_echoes.im_recon));
    if isNegPEDir
        img_all_echoes.im_recon = flip(img_all_echoes.im_recon,1);
        img_all_echoes.im_recon = flip(img_all_echoes.im_recon,2);
    end

    
    img_niiInfo.ImageSize(4:5) = img_niiInfo.ImageSize([5 4]);
    img_niiInfo.ImageSize(5) = [];
    % img_niiInfo.PixelDimensions(4:5) = img_niiInfo.PixelDimensions([5 4]);
    img_niiInfo.PixelDimensions(4) = TR;
    img_niiInfo.PixelDimensions(5) = [];
    img_niiInfo.SpaceUnits = 'Millimeter';
    img_niiInfo.TimeUnits = 'Second';
    img_niiInfo.raw.dim(1) = 4;
    img_niiInfo.raw.dim(5:6) = img_niiInfo.raw.dim([6 5]);
    % img_niiInfo.raw.pixdim(5:6) = img_niiInfo.raw.pixdim([6 5]);
    img_niiInfo.raw.pixdim(5) = TR;
    img_niiInfo.raw.xyzt_units = 10;
    
    
    nifti_fname = ['./result/im',num2str(mid(1)),'.result_',reconName,'.nii'];
    niftiwrite(img_all_echoes.im_recon,nifti_fname,img_niiInfo,'Compressed',true);


end

