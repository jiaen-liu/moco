%%% This is a demo showing how to reconstruct BOLD fMRI data. 
close all; clc; clearvars;

%% Specify the paths for the working and data directories
workDir= '~/data-sharing/moco'; %'~/myProjects/moco/';
% workDir='~/Documents/moco-recon/';
addpath(genpath(workDir));

dataDir= '~/data-sharing/moco/demo_bold/data/'; %'~/myData/moco-data/';
cd(dataDir);

% Specify the measurement number
mid_recon = 424;
mid_ref = 403;

isBOLD = 1; % convert the .mat to .nifti in results
if isBOLD
    isNegPEDir = 0; % 0--AP 1--PA
end
ifRecon = 1;

% Specify the configuration file for recon
% 1 -- motion & B0 correction
% 2 -- B0 correction
% 3 -- no correction
disp('Reading recon configuration...')
magnet = '10p5t'; % '7t'
orient = 'sag'; % 'axi'
whichReconPar = 1;
switch whichReconPar
    case 1
        reconPar = 'steParMoCoB0Co'; %motion and B0 correction
        reconid = 'moco';
        noB0Co = 0;
    case 2
        reconPar = 'steParB0Co'; %B0 correction 
        reconid = 'nomoco';
        noB0Co = 0;
    case 3
        reconPar = 'steParNoCo'; %no correction
        reconid = 'noco';
        noB0Co = 1;
    otherwise
    error('WRONG recon configuration input!');
end

reconConf = [reconPar,'_',magnet,'_',orient,'.conf'];
recon_conf_file = fullfile(workDir,'recon_conf',reconConf);


mid2 = mid_ref;
mid1 = mid_recon;

if ifRecon
    disp('******************')
    disp(['**** MID00',num2str(mid1),' ****'])
    %% Preparation 
    disp(['-> Preparing recon for MID00',num2str(mid1), '...'])
    prep_ste(mid1,'mid_pimg',mid2,'no_b0_main',noB0Co);

    %% Reconstruction
    disp(['-> Reconstructing for MID00',num2str(mid1), '...'])
    reconAMRIMoCo(dataDir,mid1,{recon_conf_file});
end

% 
if isBOLD
    disp(['-> Converting hdr to BOLD for MID00',num2str(mid1), '...'])
    volumeTR = 2.34;
    convertBOLD(mid1,volumeTR,reconid,isNegPEDir);
    i = i+1;
end




