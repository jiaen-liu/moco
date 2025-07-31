# moco
This repository contains Matlab code that can be used to reconstruct high resolution 3D T2* weighted images with navigator-based joint motion and field correction. 
The raw k-space data should be acquired using a custom-built pulse sequence that allows collection of motion robust single or multiple gradient echo data using a GRE or EPI readout. 

A more detailed description of the reconstruction and pulse sequence can be found in the following papers: 

1. Liu J, van Gelderen P, de Zwart JA, Duyn JH. Reducing motion sensitivity in 3D high-resolution T2*-weighted MRI by navigator-based motion and nonlinear magnetic field correction. Neuroimage. 2020 Feb 1;206:116332. doi: 10.1016/j.neuroimage.2019.116332. Epub 2019 Nov 2. PMID: 31689535; PMCID: PMC6981037.
2. van Gelderen P, Li X, de Zwart JA, Beck ES, Okar SV, Huang Y, Lai K, Sulam J, van Zijl PCM, Reich DS, Duyn JH, Liu J. Effect of motion, cortical orientation and spatial resolution on quantitative imaging of cortical R2* and magnetic susceptibility at 0.3 mm in-plane resolution at 7 T. Neuroimage. 2023 Apr 15;270:119992. doi: 10.1016/j.neuroimage.2023.119992. Epub 2023 Feb 27. PMID: 36858332; PMCID: PMC10278242.


## Prerequisites

### NOTE: You can only run this recon code on a Linux machine. 
To run the reconstruction code successfully, you will need to complete the following:

1. **Intel MKL Library**

   Download the tarball named `lib_mkl_tbb*.tar` under the subfolder "Intel MKL library" shared at [this google drive](https://drive.google.com/drive/folders/1cVI2BXiPV-lKmIz1KD7RiYVmy8S9kSTL?usp=drive_link).
   Please extract the files (e.g., using `tar -xvf lib_mkl_tbb_2024_2.tar`) to a subfolder named "lib_mkl_tbb" under the Intel MKL directory, that is, to `<your_path>/moco/intel_mkl/lib_mkl_tbb/`.

   If you are using a server with an AMD CPU, you can set the following environment variable before launching MATLAB to maximize Intel MKL performances:   `export MKL_DEBUG_CPU_TYPE=5`.

4. **Binary Files for Reconstruction**
    
   Download the subfolder "amriMoCo" (storing the binary files required by reconstruction) from [this google drive](https://drive.google.com/drive/folders/1cVI2BXiPV-lKmIz1KD7RiYVmy8S9kSTL?usp=drive_link), and add it to your system PATH: 
    
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
