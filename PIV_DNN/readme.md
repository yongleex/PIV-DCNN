# Details of the PIV-DCNN
 1. File descriptions
 2. Our PIV-DCNN structure
 3. The synthetic  training data
 4. Network definition
 5. Training operation
 ---

### 1. File descriptions
```
.
├── readme.md    #The name can explain everything (this file).
├── train.m      #Generating training data, and conducting optimization of each net. The final output is 'Nets.mat ', which means you  can use our results to analysis. 
├── pivdnn.m     #Our PIV-DCNN implementation with trained nets, the parameters in 'Nets.mat'
├── cnn_pivdnn_init.m    #This file gives the definition of each individual net (Net structure).
├── cnn_pivdnn.m         #This file is related to the training procedure of each net(read data, optimization)
├── genPIVImgDB.m        #You can use this file to generate large amount of training data. 
├── genPIVImgDB_F4.m     #The networs in Net4 were trained using the data generated with this file. 
├── ImgData.mat          #Two default images of 'pivdnn.m',  represents solid rotation.
└── Nets.mat             #This file stores the  parameters of 6 nets  
```


### 2. Our PIV-DCNN structure
The detail of PIV-DCNN structure is described in Fig. 3 of the manuscript.  It is a different evaluation scheme (PIV-DCNN) with multi-level regression deep convolutional neural networks. 

###  3. The synthetic  training data
$x^2$

### 4. Network definition

### 5. Training operation

### How to train the PIVnet
- Install the MatConvNet in dir "PIVnet/matconvnet/"
- run the "Train.m", it will cost more than 2 days if you employ a GPU for accelaration. You will get a file "Nets.mat", which contains the trained DCNNs.
- You can test the performance of PIVnet  using pre-trained nets if you run "PIVnet.m" directly.

## Training procedure
- generate the training data using "genPIVImgDB.m"
- train the net with "cnn_pivnet.m"
- define the net with "cnn_pivnet_init.m"
- details in "Train.m"

