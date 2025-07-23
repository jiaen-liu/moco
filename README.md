# moco
This repository contains code to correctfor rigid-body motion and B0 change effect in high resolution T2*w MRI data with simultanious volumetric navigator acqusition. You can run the script `demo_BOLD.m` to test it. 

## Prerequisites
Before running the `demo.m` script, please ensure the following prerequisites are completed. 

1. **Intel MKL Library**

    First, you should add Intel MKL library. It's too large to include on GitHub. You can locate it in `/home/naxos2-raid26/jiaen/` with the filename `lib_mkl_tbb*.tar.gz`. Please copy and unzip the files, and put them in a folder named `lib_mkl_tbb` under the Intel MKL directory. The path should be: `<your_path>/moco_master/intel_mkl/lib_mkl_tbb/`.

2. **Environment Variable for AMD CPUs** 

    If you run the code on Skopeineers that uses AMD cpu, please set the following environment variable before starting MATLAB and running the reconstruction: 
    
    `export MKL_DEBUG_CPU_TYPE=5`

3. **Binary Files for Reconstruction**
    
    During the recon process, several binary files are required, which can be located in `/home/naxos2-raid26/jiaen/bin/amriMoCo/`. Please add them to the PATH: 
    
    `export PATH=$PATH:/home/naxos2-raid26/jiaen/bin/amriMoCo/`

## Demo script
For a quick demostration, rawdata of 2 mm low-resolution 3D BOLD EPI with 10 frames are included in the Release. Please put them under the `../moco-data/` directory. Then run the `demo.m` in the directory.

First run the preparation function: 

`prep_ste(mid1, ‘mid_pimg’, mid2);`

Here, `mid1` is the MID of the main scan. `72` for the example below. `mid2` is the MID of the reference scan, `69` as shown below.

Example:

`meas_MID00072_FID04010_AMRI_epi_EPI_2mm_ipat2x3.dat`
`meas_MID00069_FID04007_AMRI_epi_SENSEcalib_4mm.dat`

Second run the reconstruction function:

`reconAMRIMoCo(‘path_to_data’, mid2, {‘recon_conf_file’});`

`path_to_data` is the path for the data directory. `recon_conf_file` is the path for the configuration file including some recon parameters. You can find it under `../moco-recon/recon_conf` of the moco package.

The recon restults will be stored in the `../moco-data/result`.


## Reference
For reference, please check the following two papers:
1. Liu J, van Gelderen P, de Zwart JA, Duyn JH. Reducing motion sensitivity in 3D high-resolution T2*-weighted MRI by navigator-based motion and nonlinear magnetic field correction. Neuroimage. 2020 Feb 1;206:116332. doi: 10.1016/j.neuroimage.2019.116332. Epub 2019 Nov 2. PMID: 31689535; PMCID: PMC6981037.
2. van Gelderen P, Li X, de Zwart JA, Beck ES, Okar SV, Huang Y, Lai K, Sulam J, van Zijl PCM, Reich DS, Duyn JH, Liu J. Effect of motion, cortical orientation and spatial resolution on quantitative imaging of cortical R2* and magnetic susceptibility at 0.3 mm in-plane resolution at 7 T. Neuroimage. 2023 Apr 15;270:119992. doi: 10.1016/j.neuroimage.2023.119992. Epub 2023 Feb 27. PMID: 36858332; PMCID: PMC10278242.
