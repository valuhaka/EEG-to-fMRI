import torch
from torch import nn
import wandb
import numpy as np
import os
import matplotlib.pyplot as plt

import geoopt

def make_complex_loss_function(
    mse_weight=0.0, corr_weight=0.0, manifold_weight=0.0, bound=1
):
    mse_loss = nn.MSELoss()
    cos_sim = nn.CosineSimilarity(dim=-1, eps=1e-08)
    spd_manifold = geoopt.SymmetricPositiveDefinite()

    def loss_func(y_hat, y_true):
        """
        y.shape [batch, roi, time]
        """
        y_hat = y_hat.to(y_true.device)  # ✅ Move y_hat to the same device as y_true
        y_true = y_true.to(y_hat.device)  # ✅ Move y_true to the same device as y_hat
        
        #print(f"y_hat shape before slicing: {y_hat.shape}")
        #print(f"y_true shape before slicing: {y_true.shape}")
        #breakpoint()
        #print(f"y_hat device: {y_hat.device}, y_true device: {y_true.device}")  # ✅ Debug print

        batch = y_hat.shape[0]
        # L1/L1 loss.
        mse = mse_loss(y_hat[..., bound:-bound], y_true[..., bound:-bound])
        #print(f"y_hat shape: {y_hat.shape}, y_true shape: {y_true.shape}, bound: {bound}")

        #print(f"mse device: {mse.device}")  # ✅ Debug print

        # Correlation
        y_hat_centre = y_hat - torch.mean(y_hat, -1, keepdim=True)
        y_true_centre = y_true - torch.mean(y_true, -1, keepdim=True)

        corrs = cos_sim(y_hat_centre, y_true_centre)
        corr = torch.mean(corrs)
        corr = torch.nan_to_num(corr, nan=0.0)

        corr_neg = -corr

        # Manifold covariance loss
        cov_matrix_hat = torch.stack([torch.cov(y_) for y_ in y_hat])
        cov_matrix_true = torch.stack([torch.cov(y_) for y_ in y_true])

        ## Manifold distance hard-coded to zero to avoid convergence errors.

        manifold_distance = torch.zeros(1, device=y_hat.device)  # ✅ Now it's on the correct device

        """ 
        man_dists = []
        for batch in range(batch):
        
           man_dist = spd_manifold.dist(x=cov_matrix_hat[batch], 
                                         y= cov_matrix_true[batch])
           man_dists.append(man_dist)
        
        
        manifold_distance = torch.mean(torch.stack(man_dists))
        manifold_distance = torch.clip(manifold_distance, 0, 100) # values might be very big
        
        manifold_distance = torch.nan_to_num(manifold_distance, nan=0.0)
        """

        # Make sure mse_weight is also a tensor on the correct device
        mse_weight_tensor = torch.tensor(mse_weight, device=y_hat.device)
        corr_weight_tensor = torch.tensor(corr_weight, device=y_hat.device)
        manifold_weight_tensor = torch.tensor(manifold_weight, device=y_hat.device)
        #print(f"  mse_weight_tensor device: {mse_weight_tensor.device}")
        #print(f"  corr_weight_tensor device: {corr_weight_tensor.device}")
        #print(f"  manifold_weight_tensor device: {manifold_weight_tensor.device}")
    
        total_loss = (
            mse_weight_tensor * mse
            + corr_weight_tensor * corr_neg
            + manifold_weight * manifold_distance
        )
        return total_loss, corr, mse, manifold_distance

    return loss_func


def make_mse_loss():
    criterion = nn.MSELoss()
    cos_metric = nn.CosineSimilarity(dim=-1, eps=1e-08)

    def loss_func(y_hat, y_batch):
        mse_loss = criterion(y_hat, y_batch)

        y_hat_centre = y_hat - torch.mean(y_hat, -1, keepdim=True)
        y_true_centre = y_batch - torch.mean(y_batch, -1, keepdim=True)

        cos_dist = torch.mean(cos_metric(y_hat_centre, y_true_centre))
        return mse_loss, cos_dist

    return loss_func


def make_mae_loss():
    criterion = nn.L1Loss()
    cos_metric = nn.CosineSimilarity(dim=-1, eps=1e-08)

    def loss_func(y_hat, y_batch):
        mae_loss = criterion(y_hat, y_batch)
        cos_dist = torch.mean(cos_metric(y_hat, y_batch))
        return mae_loss, cos_dist

    return loss_func

'''
def train_step(x_batch, y_batch, model, optimizer, loss_function, scheduler=None):
    optimizer.zero_grad()
    y_hat = model(x_batch)
    losses = loss_function(y_hat, y_batch)

    losses[0].backward()
    optimizer.step()
    if scheduler is not None:
        scheduler.step()
    return losses'''


def train_step(x_batch, y_batch, model, optimizer, loss_function, device):
    optimizer.zero_grad()
    
    x_batch = x_batch.to(device, dtype=torch.float)
    y_batch = y_batch.to(device, dtype=torch.float)  # ✅ Move target labels to the same device as model

    y_hat = model(x_batch)  # Model prediction
    y_hat = y_hat.to(device, dtype=torch.float)  # ✅ Move predictions to the correct device

    # Debugging: Print shapes
    #print(f"y_hat shape: {y_hat.shape}, y_true shape: {y_batch.shape}")
    
    losses = loss_function(y_hat, y_batch)  # ✅ Now both y_hat & y_batch are on the same device
    losses[0].backward()
    
    optimizer.step()
    return losses



def save_checkpoint_custom(state, best_model_path):
    """
    Save checkpoint based on state information. Save model.state_dict() weight of models.
    Parameters:
    state: torch dict weights
        model.state_dict()
    best_model_path: str
        path to save best model( copy from checkpoint_path)
    """

    best_check = os.path.split(best_model_path)[0]
    if not os.path.exists(best_check):
        os.makedirs(best_check)

    torch.save(state, best_model_path)

def wanb_train_regression(
    EPOCHS,
    model,
    train_loader,
    val_loader,
    loss_function,
    train_step,
    optimizer,
    device,
    raw_test_data,
    labels,
    inference_function,
    to_many,
    scheduler, #None, #None originally
    show_info=1,
    num_losses=10,
    #patience = 15,
    #delta = 0.01,
):
    """
    Train model with train_loader.
    """
    min_loss_r = float('inf')
    patience_counter = 0
    
    max_cos_val = -1
    batch_size = train_loader.batch_size

    # X_test, y_test = raw_test_data

    print(
        "Starting Training of our model",
        "\nNumber of samples", batch_size * len(train_loader),
        "\nSize of batch:", batch_size,
        "Number batches", len(train_loader),
    )

    # -----------------------------------------------------------------------#

    model = model.to(device)
    #     wandb.watch(model, loss_function, log_freq=16)

    for epoch in range(1, EPOCHS + 1):
        model = model.to(device)

        sum_losses = [0 for i in range(num_losses)]
        sum_losses_val = [0 for i in range(num_losses)]

        # model training
        model.train()
        #breakpoint()
        for counter, (x_batch, y_batch) in enumerate(train_loader):
            x_batch = x_batch.to(device, dtype=torch.float)
            y_batch = y_batch.to(device, dtype=torch.float)

            losses = train_step(x_batch, y_batch, model, optimizer, loss_function, device)

            num_losses = len(losses)

            for i in range(num_losses):
                sum_losses[i] = sum_losses[i] + losses[i].item()

            print(".", sep=" ", end="", flush=True)

        # model validation
        model.eval()
        with torch.no_grad():
            for counter_val, (x_batch_val, y_batch_val) in enumerate(val_loader):
                x_batch_val = x_batch_val.to(device, dtype=torch.float)
                y_batch_val = y_batch_val.to(device, dtype=torch.float)

                y_hat_val = model(x_batch_val)
                losses_val = loss_function(y_hat_val, y_batch_val)

                for i in range(len(losses_val)):
                    sum_losses_val[i] = sum_losses_val[i] + losses_val[i].item()

            ### add to wanb all losses.
            mean_losses = [loss / (counter + 1) for loss in sum_losses]
            mean_losses_val = [loss / (counter_val + 1) for loss in sum_losses_val]

            if scheduler is not None:
                scheduler.step()
                
            
            #log metrics to wandb
            for i in range(num_losses):
                loss_names = ["total", "r", "MSE", "d_manifold"]
                wandb.log({"train/loss_" + loss_names[i]: mean_losses[i]}, epoch)
                wandb.log({"val/loss_" + loss_names[i]: mean_losses_val[i]}, epoch)
                      
            # inference only when cosine distance imroves
            if max_cos_val < mean_losses_val[1]:
                max_cos_val = mean_losses_val[1]

                fig, fig_bars, corrs, y_hats, y_test = inference_function(model, raw_test_data, labels=labels, 
                                                          device=device, to_many=to_many)

                wandb.log({"val_viz/plot_ts_image": wandb.Image(fig)}, epoch)
                wandb.log({"val_viz/plot_corrs": wandb.Image(fig_bars)}, epoch)
                wandb.log({"val/corr_mean": np.mean(corrs)}, epoch)
                plt.close(fig)
                plt.close(fig_bars)

                # save model in that case.
                # save weights
                filename = "epoch_{}_val_corr{:.2}.pt".format(epoch, np.mean(corrs))
                filepath_name = os.path.join(wandb.run.dir, filename)
                save_checkpoint_custom(model.state_dict(), filepath_name)
                
                
            '''
            #early stopping
            val_loss_r = mean_losses_val[1]
            if val_loss_r < min_loss_r - delta: 
                min_loss_r = val_loss_r
                patience_counter = 0 #reset patience
            else:
                patience_counter += 1 #increment patiens counter if validation correlation has not improved
                
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch}.")
                print(f"val_loss_r = {val_loss_r}")
                break
            '''

        # Logging and saving
        # -------------------------------------------------------------------------------#
        if epoch % show_info == 0:
            general_out_ = "\nEpoch {} ".format(epoch)
            for i in range(len(sum_losses)):
                tmp_string = "train loss_{} : {:.3} ".format(i, mean_losses[i])
                tmp_string_val = "val loss_{} : {:.3} ".format(i, mean_losses_val[i])
                general_out_ = general_out_ + tmp_string + tmp_string_val
            print(general_out_)

        #val_loss = mean_losses_val[0]
        #val_acc = mean_losses_val[1]
    
    return model

