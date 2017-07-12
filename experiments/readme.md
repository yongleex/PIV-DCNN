# The arrangement of experimental code

- **exp1.m**: A simple demo for PIV-DCNN.

 1. **Inputs:**   Any two particle images ( pair)
 2. **Methods:**  FFTCC, WIDIM, PIV-DCNN
 3. **Evaluation Criteria:** Visualized vector field, Magnitude contour , (u+v) histogram,  (turbulent energy spectrum).

|
|:---------------|:---------------|
|**Inputs:**| Any two particle images ( pair)| |
|**Methods:**|FFTCC, WIDIM, PIV-DCNN|
|**Evaluation Criteria:**|Visualized vector field, Magnitude contour , (u+v) histogram,  (turbulent energy spectrum)|


- **exp2.m**: The performance investigation (RMS error) of PIV-DCNN with respect to different **particle diameters**. 

	|
	|:---------------|:---------------|:-------------|
	|**Inputs:**| synthetic particle images with uniform flow |  several particle diameters|
	|**Methods:**|PIV-DCNN|
	|**Evaluation Criteria:**| RMSE curve, (Mean bias curve)|

- **exp3.m**: The performance investigation (RMS error) of PIV-DCNN with respect to different **particle concentration** (unit: particles per pixel,ppp).

	|
	|:---------------|:---------------|:-------------|
	|**Inputs:**| synthetic particle images with uniform flow |  several concentrations|
	|**Methods:**|PIV-DCNN|
	|**Evaluation Criteria:**| RMSE curve, (Mean bias curve)|

- **exp4.m**: The performance investigation (RMS error) of PIV-DCNN with respect to different **noise level**. 

	|
	|:---------------|:---------------|:-------------|
	|**Inputs:**| synthetic particle images with uniform flow |  several noise levels|
	|**Methods:**|PIV-DCNN|
	|**Evaluation Criteria:**| RMSE curve, (Mean bias curve)|

- **exp5.m**: The performance investigation (RMS error) of **PIV-DCNN** in comparison with **FFT-CC** and **WIDIM** methods. 

	|
	|:---------------|:---------------|
	|**Inputs:**| synthetic particle images with uniform flow |
	|**Methods:**|FFTCC, WIDIM, PIV-DCNN|
	|**Evaluation Criteria:**| RMSE curve, (Mean bias curve)|

- **exp6.m**: The Modulation Transfer Function of PIV-DCNN, i.e, **spatial resolution test**. 

	|
	|:---------------|:---------------|
	|**Inputs:**| synthetic particle images with sinusoidal flow |  
	|**Methods:**|FFTCC, WIDIM, PIV-DCNN|
	|**Evaluation Criteria:**|MC response curve, (RMSE curve)|

- **exp7.m**: The filters visualization.  

 



 

