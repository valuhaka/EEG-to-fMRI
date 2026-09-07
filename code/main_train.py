import sys, os
import numpy as np
import pandas as pd
from datetime import datetime
import wandb
#import matplotlib.pyplot as plt

import torch
import torch.optim as optim
#from torch.utils.data import Subset
#from pytorch_model_summary import summary

sys.path.insert(1, os.path.realpath(os.path.pardir))
os.chdir("/data/p_02950/restore_20250928-220205/scripts/BEIRA/3_master/") #

from utils import preproc
from utils import torch_dataset
from utils import train_utils
from utils import inference
from utils.models_arch import autoencoder_new_ArturNH

##--------------------------------------------------------------------------##
base_directory = "/.../" #insert path to where you stored the repository
perm_base_directory = os.path.join(base_directory, "data", "DL_model")

permutation_folders = [folder for folder in os.listdir(perm_base_directory) 
                      if os.path.isdir(os.path.join(perm_base_directory, folder))]
permutation_folders.sort()  # Sort for consistent ordering
print(f"Found permutation folders: {permutation_folders}")

roi_options = dict(
        VS = ["VS"],
        Kovalev8=["Left Caudate", "Left Putamen", "Left Pallidum", "Left Accumbens", "Right Caudate", "Right Putamen", "Right Pallidum", "Right Accumbens"]) 

rois = 'VS'
sub_n = 'train_set.npz' 
project_name = "FINAL"   #name of project in wandb
device = "cuda" if torch.cuda.is_available() else "cpu"


def print_input_hook(module, input, output):
    if not hasattr(print_input_hook, "printed"):  # Check if we've printed before
        print("Input to the first layer:", input[0].shape)
        print_input_hook.printed = True  # Set the attribute to prevent printing again


def fit(sub_n, rois, perm_folder):
    # Set experiment name to include the permutation folder
    exp_name = f"{perm_folder}_"
    
    seed = 42#np.random.randint(1, 1000)
    torch.manual_seed(seed) #set fixed seed for reproducibility for HP tuning
    
    config = dict(
        subject=sub_n,
        permutation=perm_folder,  # Add permutation info to config
        wandb_path=os.path.join(base_directory, "wandb_logs/"),
        n_fits=1,  
        n_channels = 62,                       
        chosen_rois = roi_options[rois],        
        n_runs=4,
        eeg_final_fq=100,
        data_limit=5,   
        freqs=[-1],    
       
        bold_delay=6,
        to_many=True,
        random_subsample=True,
        sample_per_epoch=512, 
        WINDOW_SIZE=16384,  #256,#512,#1024,#2048,#4096,#8192 #16384,
        
        lr=1e-6,
        weight_decay=5e-2,
        batch_size=32,
        
        mse_weight=0.9,
        corr_weight=0.1,
        manifold_weight=0,
)
       
    if config["n_channels"] == 29:
        config["eeg_channels"] = ["Fp1", "Fp2",                 # in Kovelev "Fz" -> average from Fp1 & Fp2 to form Fz?
                              "F7", "F3", "Fz", "F4", "F8", "FC5", "FC6", "T7", "C3", "Cz", "C4", "T8", "TP9", "CP5", "CP6", "TP10", "P7", "P3", "Pz", "P4", "P8",
                              "PO3", "POz", "PO4", "PO9", "Oz", "PO10"]         #no Iz in Andreou data - so only 29 electrodes 
        
    elif config["n_channels"] == 62:
        config["eeg_channels"] =   ["Fp1", "Fp2", "F3", "F4", "C3", "C4",
                                "P3", "P4", "O1", "O2", "F7", "F8", "T7", "T8", "P7", "P8", "Fz", "Cz", "Pz", "Oz", "FC1", "FC2", "CP1", "CP2", "FC5", "FC6", "CP5", "CP6", "TP9", "TP10",
                                "F1", "F2", "C1", "C2", "P1", "P2", "AF3", "AF4", "FC3", "FC4", "CP3", "CP4", "PO3", "PO4", "F5", "F6", "C5", "C6", "P5", "P6", "AF7", "AF8", 
                                "FT7", "FT8", "TP7", "TP8", "FT9", "FT10", "PO9", "PO10", "CPz", "POz"]
       
    hp_autoencoder = dict(
        n_electrodes=len(config["eeg_channels"]),
        n_channels_out=len(config["chosen_rois"]), #why channels out? I think this should be roi_out
        n_freqs=len(config["freqs"]),
        channels=[128, 128, 128, 128], #[8, 16, 32, 64], #[128, 128, 128, 128] originally
        kernel_sizes=[3,3,3], # 3, 5, 7 in the previous model
        strides=[8,8,8],
        dilation=[1,1,1],
        decoder_reduce=4,
        dropout_rate = 0.4, #0.3 originally
        hidden_channels=16, #originally 16
    )

    config = {**hp_autoencoder, **config}  
    
    # Update train_directory to use the specific permutation folder
    train_directory = os.path.join(perm_base_directory, perm_folder)
    files = os.listdir(train_directory)                                #make list of files in the preproc folder  
    subfiles = [file for file in files if file.endswith('.npz')]
    
    #load train data 
    sub_path = os.path.join(train_directory, sub_n)  #path to where the training data is
    train_data = np.load(sub_path, allow_pickle=True)
    train_eeg = train_data['eeg']
    
    #extract the relevant channels
    eeg_channel_names = train_data['eeg_channel_names']  # List of all EEG channel names
    channel_inds = [eeg_channel_names.tolist().index(ch) for ch in config["eeg_channels"] if ch in eeg_channel_names] # Find indices of the channels in eeg_channel_names that match config_channels
    train_eeg = train_eeg[channel_inds, :]         # Filter the EEG data to retain only the specified channels
    
    train_fmri = train_data['fmri'][1:2, :]
          
    #load test data  
    sub_path = os.path.join(train_directory, "test_set.npz")  #path to where the test data is
    test_data = np.load(sub_path, allow_pickle=True)
    test_eeg = test_data['eeg']
    eeg_channel_names = test_data['eeg_channel_names']  # List of all EEG channel names
    channel_inds = [eeg_channel_names.tolist().index(ch) for ch in config["eeg_channels"] if ch in eeg_channel_names] # Find indices of the channels in eeg_channel_names that match config_channels
    test_eeg = test_eeg[channel_inds, :]         # Filter the EEG data to retain only the specified channels
    
    test_fmri = test_data['fmri'][1:2, :]
    
    #NORMALIZATION
    # Normalize EEG data (per channel, across time)
    train_eeg_mean = np.mean(train_eeg, axis=1, keepdims=True)
    train_eeg_std = np.std(train_eeg, axis=1, keepdims=True)
    train_eeg = (train_eeg - train_eeg_mean) / (train_eeg_std + 1e-8)
    
    # Normalize fMRI data (per ROI, across time)
    train_fmri_mean = np.mean(train_fmri, axis=1, keepdims=True)
    train_fmri_std = np.std(train_fmri, axis=1, keepdims=True)
    train_fmri = (train_fmri - train_fmri_mean) / (train_fmri_std + 1e-8)
    
    # Apply same normalization parameters to test data
    test_eeg = (test_eeg - train_eeg_mean) / (train_eeg_std + 1e-8)
    test_fmri = (test_fmri - train_fmri_mean) / (train_fmri_std + 1e-8)
   
    train_dataset_prep = (train_eeg, train_fmri)
    test_dataset_prep = (test_eeg, test_fmri)


    # Apply time delay.
    train_dataset_prep = preproc.bold_time_delay_align(
        train_dataset_prep, config["eeg_final_fq"], config["bold_delay"]
    )

    test_dataset_prep = preproc.bold_time_delay_align(
        test_dataset_prep, config["eeg_final_fq"], config["bold_delay"]
    )

    print("Size of train dataset:", train_dataset_prep[0].shape, train_dataset_prep[1].shape)
    print("Size of test dataset:", test_dataset_prep[0].shape, test_dataset_prep[1].shape)

    # Torch dataset creation.
    torch_dataset_train = torch_dataset.CreateDataset_eeg_fmri(
        train_dataset_prep,
        random_sample=config["random_subsample"],
        sample_per_epoch=config["sample_per_epoch"],
        to_many=config["to_many"],
        window_size=config["WINDOW_SIZE"],
    )

    torch_dataset_test = torch_dataset.CreateDataset_eeg_fmri(
        test_dataset_prep,
        random_sample=False,
        sample_per_epoch=None,
        to_many=config["to_many"],
        window_size=config["WINDOW_SIZE"],
    )
    print("Size of test dataset:", len(torch_dataset_test))

    params_train = {"batch_size": config["batch_size"], "shuffle": False, "num_workers": min(8, os.cpu_count())}
    params_val = {"batch_size": config["batch_size"], "shuffle": False, "num_workers": min(4, os.cpu_count())}

    # init dataloaders for training
    train_loader = torch.utils.data.DataLoader(torch_dataset_train, **params_train)
    val_loader = torch.utils.data.DataLoader(torch_dataset_test, **params_val)


    ##-----------------------------------------------------------------------##
    ## Fit it.
    time_log = pd.Series()
    time_logs = pd.DataFrame()
       
   
    print("Training the model for", config["subject"], "on", perm_folder)

    for ind in range(config["n_fits"]):
        start_time = datetime.now()
        time_log["start time"] = str(start_time.time())

        print("\nRun", ind, "initiated at:", start_time.time(), "\n")
               
        # Model initialization
        model = autoencoder_new_ArturNH.AutoEncoder1D_Artur_MultiHead(hp_autoencoder)
        if torch.cuda.device_count() > 1:
            print(f"Using {torch.cuda.device_count()} GPUs!")
            model = torch.nn.DataParallel(model)  # Wrap model in DataParallel
            
        model.to(device)
        #print(model)
        
        #model.load_state_dict(state_dict, strict=True)      # Load the weights 

        loss_func = train_utils.make_complex_loss_function(
            mse_weight=config["mse_weight"],
            corr_weight=config["corr_weight"],
            manifold_weight=config["manifold_weight"],
            bound=1,
        )

        train_step = train_utils.train_step
        
        optimizer = optim.AdamW(model.parameters(), 
                                lr=config["lr"], 
                                weight_decay=config["weight_decay"])
        
        
        from torch.optim.lr_scheduler import CosineAnnealingLR
        scheduler = CosineAnnealingLR(optimizer, T_max=500)

        
        parameters = {
            "EPOCHS": 1500,
            "model": model,
            "train_loader": train_loader,
            "val_loader": val_loader,
            "loss_function": loss_func,
            "train_step": train_step,
            "optimizer": optimizer,
            "device": device, 
            "raw_test_data": test_dataset_prep,
            "show_info": 1,
            "num_losses": 5,
            "labels": config["chosen_rois"],
            "inference_function": inference.model_inference_function,
            "to_many": config["to_many"],
            "scheduler": scheduler
        }

        path_to_save_wandb = config["wandb_path"]
        
        wandb.login(key='insert_key') 
        
        with wandb.init(
            project= project_name, 
            config=config, 
            save_code=True, 
            dir=path_to_save_wandb
        ):
            wandb.define_metric("val/corr_mean", summary="max")

            if ind == 0:
                wandb.run.name = exp_name + wandb.run.name #config['subject'][0:6]
                            
            model = train_utils.wanb_train_regression(**parameters)
            model.to(device)
	    
            # Save the weights with permutation info.
            weights_name = f"weights_from_{wandb.run.name}_run_{ind}.pth"

            torch.save(
                model.state_dict(), os.path.join(path_to_save_wandb, "weights", weights_name)
            )
            
            ##---------------------------------------------------------------##
            ## Table and print the time spent in the logs.

            end_time = datetime.now()
            time_log["end_time"] = str(end_time.time())

            time_diff = end_time - start_time
            time_diff_per_epoch = time_diff / parameters["EPOCHS"]
            time_log["duration"] = str(time_diff)
            time_log["d / epoch"] = str(time_diff_per_epoch)

            print("\nRun", ind, "finished at:", end_time.time(), "\n")

            time_logs.insert(ind, ind, time_log)
            print(time_logs)

            print(
                "\n\nMSE weight:",
                config["mse_weight"],
                "\ncorr weight:",
                config["corr_weight"],
                "\nEpochs:",
                parameters["EPOCHS"],
            )
    #hook.remove()
    wandb.finish()    

##---------------------------------------------------------------------------##
# Main loop to train on multiple permutation folders
for perm_folder in permutation_folders:
    print(f"\n{'='*60}")
    print(f"Starting training for {perm_folder}")
    print(f"{'='*60}\n")
    
    fit(sub_n, rois, perm_folder)
    print(f"Successfully completed training for {perm_folder}")

print("\nAll permutation trainings completed!")
