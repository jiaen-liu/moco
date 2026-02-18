%%% This is a demo showing how to reconstruct multi-echo GRE data. 
close all; clc; clearvars;

%% Specify the paths for the working and data directories
workDir='~/data-sharing/moco'; % '~/myProjects/moco/';
addpath(genpath(workDir));

dataDir = '~/data-sharing/moco/demo_megre/data/'; %'~/myData/megre-data';
cd(dataDir);

%% Specify the measurement number
mid_recon = 547;
mid_ref = 552;

% Specify the configuration file for recon
% 1 -- motion & B0 correction
% 2 -- B0 correction
% 3 -- no correction
% 4 -- motion correction 
disp('Reading recon configuration...')
magnet = '7t'; % '10p5t'
orient = 'axi'; % 'sag'
whichReconPar = 1;
switch whichReconPar
    case 1
        reconPar = 'steParMoCoB0Co'; % motion and B0 correction
        noB0Co = 0;
    case 2
        reconPar = 'steParB0Co'; % B0 correction 
        noB0Co = 0;
    case 3
        reconPar = 'steParNoCo'; % no correction
        noB0Co = 1;
    case 4
        reconPar = 'steParMoCo'; % motion correction 
        noB0Co = 1;
    otherwise
    error('WRONG recon configuration input!');
end

reconConf = [reconPar,'_',magnet,'_',orient,'.conf'];
recon_conf_file = fullfile(workDir,'recon_conf',reconConf);


mid2 = mid_ref;
mid1 = mid_recon;

disp('******************')
disp(['**** MID00',num2str(mid1),' ****'])

%% Preparation 
disp(['-> Preparing recon for MID00',num2str(mid1), '...'])
prep_ste(mid1,'mid_pimg',mid2,'no_b0_main',noB0Co);

%% Reconstruction
disp(['-> Reconstructing for MID00',num2str(mid1), '...'])
reconAMRIMoCo('.',mid1,{recon_conf_file});






