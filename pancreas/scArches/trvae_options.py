#!/PHShome/ik936/anaconda3/envs/gpu/bin/python

## SETTINGS
# trvae_epochs = 1
# surgery_epochs = 1
# nepochs_alpha = 1
trvae_epochs = 500
surgery_epochs = 500
nepochs_alpha = 200

early_stopping_kwargs = {
    "early_stopping_metric": "val_unweighted_loss",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
