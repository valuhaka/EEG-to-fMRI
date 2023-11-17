# EEG-to-fMRI
Predicting reward-related BOLD activation from EEG data.


# EEG

1. Preprocessing

- interpolation
- ICA
- **no** cutting


# fMRI

1. Preprocessing

- head movement
- selecting striatal ROI

2. Fitting

- hemodynamic response function


## Modeling

- autoencoder

1. Finding latent factors from EEG

2. Generating data

3. Comparing generated data to fMRI

