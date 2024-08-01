# Compiling and Running mSWEEP-GPU on HPC environments
## Turso HPC (University of Helsinki)

Has NVIDIA A100 GPUs with 80GB memory, so CUDA LibTorch required (download and unzip somewhere)

### Before compiling
Load modules before compiling:
```
module load CUDA GCC/12.2.0 cmake
```
### Before running
When running in a batch script, make sure to **request a GPU** and have **`CUDA` and `GCC` modules loaded in**.

## LUMI

Has AMD MI250x GPUs with 128GB memory, so ROCm LibTorch required (download and unzip somewhere)

### Before compiling
Load modules before compiling:
```
module load LUMI/23.09 partition/G gcc rocm/5.4.6
export PYTORCH_ROCM_ARCH=gfx90a
```

Even though there is ROCm 5.6.1 on LUMI, it doesn't seem to work properly, so 5.4.6 is the highest version which works (LibTorch version 2.0.1 for ROCm 5.4.2 works with this).

### Before running
When running in a batch script, make sure to **request a GPU** and have **`LUMI/23.09`, `partition/G`, `gcc`, and `rocm/5.4.6` modules loaded in**.
