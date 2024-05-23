# Adaptive Parameter Learning Algorithm

## Software configuration requirements
MATLAB R2022b

## Code Structure
- `Main/Main.m`: Master file
- `GAMP/Complex_OTFS_GAMP*.m`: Message passing algorithm
- `Operator/System_Complex_OTFS_Channel_Estimation.m`: Generative system model
- `Operator/Estimate_EM_ELBO_enhance.m`: Proposed algorithm
- `Operator/Estimate_EM.m`: SR-TSP15 algorithm
- `Operator/Estimate_AWGN_noise_EM.m`: VS-TSP13 algorithm
- `Input/In_Complex_Constellation_Estimation.m`: Constellation point correlation
- `Output/Out_Complex_Quantization_Estimation_Rewrite.m`: Coarsely quantized ADC correlation

## Citation
```tex
@article
```