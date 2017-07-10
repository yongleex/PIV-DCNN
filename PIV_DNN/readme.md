# How to train the PIVnet
- Install the MatConvNet in dir "PIVnet/matconvnet/"
- run the "Train.m", it will cost more than 2 days if you employ a GPU for accelaration. You will get a file "Nets.mat", which contains the trained DCNNs.
- You can test the performance of PIVnet  using pre-trained nets if you run "PIVnet.m" directly.

## Training procedure
- generate the training data using "genPIVImgDB.m"
- train the net with "cnn_pivnet.m"
- define the net with "cnn_pivnet_init.m"
- details in "Train.m"

