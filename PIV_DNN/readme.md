# Details of the PIV-DCNN
 1. File descriptions
 2. Our PIV-DCNN structure
 3. Network definition
 4. The synthetic  training data
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

### 3. Network definition
Each net is a regression net, whose structure is similar to the [LeNet -5](http://yann.lecun.com/exdb/lenet/) for feature extraction, and  the output of our net is continuous velocity component instead of classification label. Details of the Net can be found in Fig. 4 of our manuscript.

### 4. The synthetic  training data

#### 4.1 Synthetic underlying flow (vector field) in a small patch
The key problem of synthetic data is the vector field. We combine random weighted  algebra functions    to approximate the flow vector distribution in a small patch.

- $u(x,y)= \sum p_i \cdot \textbf{e}_i (x,y)$
- $v(x,y)=\sum q_i   \cdot \textbf{g}_i(x,y)$

where  p and q are random weights, and e(x,y),g(x,y) denote  the basic algebra functions.  In our implementation, we use a normalized coordinates x, y in [-0.5,0.5] range.  The specification of our implementation is below, and you can check it with our  Matlab Code.

##### 4.1.1. The  u component

|order i| e(x,y)  |  range of p | others|
|:------:|:--------|:-------------|:---------|
|1          | 1           |  [-Vmax, Vmax]              | Vector component  |
|2          | y           | [-10,10]       | -   |
|3          | y^2      | [-10,10]       | -   |
|4          | y^3      | [-10,10]       | -   |
|5          | sin(ay+b)  -sin(b)  |  [-3.2,3.2]      | a in [0,2pi], b in [0,2pi]  |
|6          | x           | [-10,10]  | &1  |
|7          | x^2      |  [-10,10] | -   |
|8          | x*y          | [-10,10]  | -   |
|9          | x^3           | [-10,10]   | -   |
|10        | x*y^2           | [-10,10]   | -   |
|11        | y*x^2           |  [-10,10]    | -   |
|12        | sin(cx+d)\*cos(cy+d)  - sin(d)\*cos(d)          | [-0.25,0.25]       | c in [0,2pi], d in [0,2pi]    |

 

##### 4.1.2. The  v component
|order i| g(x,y)  |  range of q | others|
|:------:|:--------|:--------------|:--------|
|1          | 1           | [-Vmax, Vmax]          |  vector component   |
|2          | x           | [-10,10]        | -   |
|3          | x^2       | [-10,10]       | -   |
|4          | x^3       | [-10,10]       | -   |
|5          | sin(rx+s)  - sin(s)         |  [-3.2,3.2]               | r in [0,2pi], s in [0,2pi]    |
|6          | y           | [-10,10]        | related to p6    |
|7          | x*y           |  [-10,10]   | related to  p7  |
|8          | y^2          | [-10,10]    |  related to p8   |
|9          | y*x^2           |  [-10,10]   | related to p9 |
|10        | x*y^2           | [-10,10]    | related to p10  |
|11        | y^3           | [-10,10]        | related to p11  |
|12        | cos(cx+d)\*sin(cy+d)  - cos(d)\*sin(d)         | [-0.25,0.25]      |   related to p12   |

##### 4.1.3.  Our matlab implementation
```Matlab
    u = @(x,y) P(1)+ P(2)*y+ P(3)*y.^2 + P(4)*y.^3  +  P(5)*sin(P(6)*y  +P(7))-P(5)*sin(P(7)) + P(15)*x + 0.5*P(16)*x.^2 +     P(17)*x.*y + ...
        P(18)*x.*x.*x/3 + P(19)*x.*y.*y + P(20)*x.*x.*y/2 + P(21).*sin(P(22)*x+P(23)).*cos(P(22)*y+P(23)) - P(21).*sin(P(23)).*cos(P(23));
    v = @(x,y) P(8)+ P(9)*x+P(10)*x.^2 + P(11)*x.^3 + P(12)*sin(P(13)*x+P(14))-P(12)*sin(P(14)) - P(15)*y -     P(16)*x.*y - 0.5*P(17)*y.^2 - ...
        P(18)*x.*x.*y - P(19)*y.^3/3  - P(20)*x.*y.*y/2 - P(21).*cos(P(22)*x+P(23)).*sin(P(22)*y+P(23))+ P(21).*cos(P(23)).*sin(P(23));
```
The synthetic vector fields can be found in our manuscript.


#### 4.2 Two way to generate training images with the  synthetic flow: the particles model and image warping operation 

### 5. Training operation


