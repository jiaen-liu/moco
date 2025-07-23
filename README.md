# moco
This repository contains code to correctfor rigid-body motion and B0 change effect in high resolution T2*w MRI data with simultanious volumetric navigator acqusition. You can run the script `demo_BOLD.m` to test it. 


## Prerequisites
Before running the `demo_BOLD.m` script, please ensure the following prerequisites are completed. 

Some required files and demo data can be downloaded from: 
https://drive.google.com/drive/folders/1cVI2BXiPV-lKmIz1KD7RiYVmy8S9kSTL?usp=drive_link

1. **Intel MKL Library**

    First, download the Intel MKL library file named `lib_mkl_tbb*.tar.gz` from the link above. You can access it through the Google Drive link above. Please unzip the files and place them in a folder named `lib_mkl_tbb` under the Intel MKL directory. The path should be: `<your_path>/moco/intel_mkl/lib_mkl_tbb/`.

2. **Environment Variable for AMD CPUs** 

    If you are using a server with an AMD CPU, set the following environment variable before launching MATLAB and running the reconstruction: 
    
    `export MKL_DEBUG_CPU_TYPE=5`

3. **Binary Files for Reconstruction**
    
    Several binary files are required during reconstruction. Please add them to your system PATH: 
    
    `export PATH=$PATH:<your_path>/amriMoCo/`

## Demo script
For a quick demostration, download the 2 mm low-resolution 3D BOLD EPI dataset with 10 frames from the link above. Place the data under the `../moco-data/` directory, then run the `demo_BOLD.m` script in the directory.

First run the preparation function: 

`prep_ste(mid1, ‘mid_pimg’, mid2);`

Here, `mid1` is the MID of the main scan. `72` for the example below. `mid2` is the MID of the reference scan, `69` as shown below.

Example:

`meas_MID00072_FID04010_AMRI_epi_EPI_2mm_ipat2x3.dat`
`meas_MID00069_FID04007_AMRI_epi_SENSEcalib_4mm.dat`

Second, run the reconstruction function:

`reconAMRIMoCo(‘path_to_data’, mid2, {‘recon_conf_file’});`

- `path_to_data`: Path to the raw data directory.
- `recon_conf_file`: Path to the reconstruction configuration file, which includes relevant parameters. You can find it under `../moco-recon/recon_conf` in this repository.

The recon restults will be saved under:

 `../moco-data/result`.


## Reference
For reference, please check the following two papers:
1. Liu J, van Gelderen P, de Zwart JA, Duyn JH. Reducing motion sensitivity in 3D high-resolution T2*-weighted MRI by navigator-based motion and nonlinear magnetic field correction. Neuroimage. 2020 Feb 1;206:116332. doi: 10.1016/j.neuroimage.2019.116332. Epub 2019 Nov 2. PMID: 31689535; PMCID: PMC6981037.
2. van Gelderen P, Li X, de Zwart JA, Beck ES, Okar SV, Huang Y, Lai K, Sulam J, van Zijl PCM, Reich DS, Duyn JH, Liu J. Effect of motion, cortical orientation and spatial resolution on quantitative imaging of cortical R2* and magnetic susceptibility at 0.3 mm in-plane resolution at 7 T. Neuroimage. 2023 Apr 15;270:119992. doi: 10.1016/j.neuroimage.2023.119992. Epub 2023 Feb 27. PMID: 36858332; PMCID: PMC10278242.
