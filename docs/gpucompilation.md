# Compiling and Running mSWEEP-GPU on HPC environments
## Turso HPC (University of Helsinki)

Has NVIDIA A100 GPUs with 80GB memory, so CUDA LibTorch required (download and unzip somewhere)

### Before compiling
Load modules before compiling:
```
module load CUDA GCC/12.2.0 cmake
```

Use GCC/12.2.0 instead of GCC which loads version 13.2.0, since there are some incompatibility issues with the versions of CUDA on the system and the new version of GCC that cause compilation not to work.

### Before running
When running in a batch script, make sure to **request a GPU** and have **`CUDA` and `GCC` modules loaded in**.

## LUMI

Has AMD MI250x GPUs with 128GB memory, so ROCm LibTorch required (download and unzip somewhere)

### Before compiling
Load modules and export environment variable before compiling:
```
module load LUMI/23.09 partition/G gcc rocm/5.4.6
export PYTORCH_ROCM_ARCH=gfx90a
```

Even though there is ROCm 5.6.1 on LUMI, it causes issues while compiling, so 5.4.6 is the highest version which works (LibTorch version 2.0.1 for ROCm 5.4.2 works with this).

### Before running
When running in a batch script, make sure to **request a GPU** and have **`LUMI/23.09`, `partition/G`, `gcc`, and `rocm/5.4.6` modules loaded in**.
