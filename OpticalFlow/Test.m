  load opticalFlowTest;
  [Vx,Vy]=opticalFlow(I1,I2,'smooth',1,'radius',10,'type','LK');
  figure(1); im(I1); figure(2); im(I2);
  figure(3); im([Vx Vy]); colormap jet;