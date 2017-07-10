%Testing calculation speeds for PIVlab. Feel free to share your results and your configuration on the PIVlab forums (http://pivlab.blogspot.de/p/forum.html)
clear all
close all
clc
disp('This script performs the four most time consuming calculations in PIVlab')
disp('and measures the time to complete them.')
disp('...testing calculation speed, please wait...')
%%
disp('...testing DCC...')
A=rand(64,64)*255;
B=rand(128,128)*255;
tic
for i=1:250
    %C= conv2(B,rot90(conj(A),2),'valid'); %full line of code
    Atemp=conj(A);
    Atemp=rot90(Atemp,2);
    C=conv2(B,Atemp,'valid');
end
dccspeed=toc/i*1000;
%%
disp('...testing DFT...')
A=round(rand(64,64)*255);
B=round(rand(64,64)*255);
tic
for i=1:10000
    %C =fftshift(real(ifft2(conj(fft2(A)).*fft2(B)))); %full line of code
    Atemp=fft2(A);
    Atemp=conj(Atemp);
    Btemp=fft2(B);
    Atemp=(Atemp.*Btemp);
    Atemp=ifft2(Atemp);
    Atemp=real(Atemp);
    C=fftshift(Atemp);
end
dftspeed=toc/i*1000;
%%
disp('...testing linear window deformation...')
A=679:743;
B=(71:135)';
C=round(rand(65,65)*255);
D=repmat(680:743,64,1)+rand(64,64);
E=repmat((72:135)',1,64)+rand(64,64);
tic
for i = 1:6000
    F=interp2(A,B,double(C),D,E,'*linear');
end
linspeed=toc/i*1000;
%%
disp('...testing spline window deformation...')
tic
for i = 1:1600
    F=interp2(A,B,double(C),D,E,'*spline');
end
splspeed=toc/i*1000;
disp('...finished')
%%
disp('----------------------------------------')
disp('Your results (time per operation):')
disp(['DCC calculation speed:      ' num2str(dccspeed) ' ms'])
disp(['DFT calculation speed:      ' num2str(dftspeed) ' ms'])
disp(['Linear interpolation speed: ' num2str(linspeed) ' ms'])
disp(['Spline interpolation speed: ' num2str(splspeed) ' ms'])
%% Williams results
%{
Surface Pro 3
Intel Core i5-4300U @ 1.9GHz 2.5GHz
8 GB RAM
64 bit Win 8.1 Pro
MATLAB R2014b 64 bit

DCC calculation speed:      13.8081 ms
DFT calculation speed:      0.30431 ms
Linear interpolation speed: 0.56215 ms
Spline interpolation speed: 1.9602 ms



PC
Intel Core i7-2600K @ 3.4GHz
16 GB RAM
64 bit Win 7 Pro
MATLAB R2014b 64 bit

DCC calculation speed:      6.0144 ms
DFT calculation speed:      0.28481 ms
Linear interpolation speed: 0.53755 ms
Spline interpolation speed: 1.1093 ms
%}