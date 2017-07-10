
J = (imread('A_00001_a.tif'));
J = single(J(158:1116,1240:2468))./65535;
J = imadjust(J);
imwrite(uint8(J*255),'A9.bmp');




J = (imread('B_00001.tif'));
J = single(J(300:748,81:1189))./65535;
J = imadjust(J);
imwrite(uint8(J*255),'A10.bmp');

I = imread('A1.tif')
imwrite(uint8(I),'A1.bmp');

I = rgb2gray(imread('A11.jpg'));
imwrite(uint8(I(1:8:end,1:8:end)),'A11.bmp');
imtool(I(1:8:end,1:8:end))



