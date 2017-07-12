# The arrangement of experimental code

- **exp1.m**: A simple demo for PIV-DCNN.

	 1. **Inputs:**   Any two particle images ( pair)
	 2. **Methods:**  FFTCC, WIDIM, PIV-DCNN
	 3. **Evaluation Criteria:** Visualized vector field, Magnitude contour , (u+v) histogram,  (turbulent energy spectrum).

----
- **exp2.m**: The Monte-Carlo    investigation (RMS error) of PIV-DCNN  performance with respect to different **particle diameters**. 

	 1. **Inputs:** synthetic particle images with uniform flow [several particle diameters]
	 2. **Methods:** PIV-DCNN
	 3. **Evaluation Criteria:** RMSE curve, (Mean bias curve)

----

- **exp3.m**: The Monte-Carlo    investigation (RMS error) of PIV-DCNN  performance with respect to different **particle concentration** (unit: particles per pixel,ppp).

	 1. **Inputs:** synthetic particle images with uniform flow  [several concentrations]
	 2. **Methods:** PIV-DCNN
	 3. **Evaluation Criteria:** RMSE curve, (Mean bias curve).

----

- **exp4.m**: The Monte-Carlo    investigation (RMS error) of PIV-DCNN  performance with respect to different **noise level**. 

	 1. **Inputs:** synthetic particle images with uniform flow [several noise levels]
	 2. **Methods:** PIV-DCNN
	 3. **Evaluation Criteria:** RMSE curve, (Mean bias curve)

----

- **exp5.m**: The Monte-Carlo    investigation (RMSE performance) of **PIV-DCNN** in comparison with **FFT-CC** and **WIDIM** methods. 

	1. **Inputs:** synthetic particle images with uniform flow
	2. **Methods:** FFTCC, WIDIM, PIV-DCNN
	3. **Evaluation Criteria:** RMSE curve, (Mean bias curve)

----

- **exp6.m**: The Modulation Transfer Function of PIV-DCNN, i.e, **spatial resolution test**.  Monte-Carlo  Simulations.

	 1. **Inputs:** synthetic particle images with sinusoidal flow 
	 2. **Methods:** FFTCC, WIDIM, PIV-DCNN
	3. **Evaluation Criteria:** MC response curve, (RMSE curve)

----

- **exp7.m**: The filters visualization.  

 



 

