%% Experiments 1
%- Apply PIVdnn method to evaluating the laboratory Particle images with state-of-art PIV methods
%- Yong Lee (leeyong@hust.edu.cn)
%- 2016-12-18
%- Modified @ 2017-05-19

function exp1()
close all; clear
%- Configuration
opts.show = true;
step_size = 16;

%- Add methods pathes
addpath('../PIV_DNN/','../PIVCC/','../pivlab/','../OpticalFlow/');

%- Read the images
oldPath = pwd;
cd('../data/TestImages');
[FileName,PathName,~]=uigetfile({'*.tif';'*.bmp';'*.jpg';'*.png'},'Please open 2 images!','MultiSelect','on');
FileName = sortrows(FileName);   
cd(oldPath);

I1 = imread([PathName,FileName{1}]); I2 = imread([PathName,FileName{2}]);

if isa(I1,'uint16') && isa(I2,'uint16'),  I1 = uint8(single(imadjust(I1))/257); I2 = uint8(single(imadjust(I2))/257);end % 16bit ->8 bit image
if numel(size(I1))>2, I1 = rgb2gray(I1); I2 = rgb2gray(I2); end % rgb images -> gray scale images

%- Process the images with WIndow Deformation Iterative Multigrid (WIDIM) PIV
 opt1.x_start = 5;opt1.y_start = 5; opt1.x_step = step_size;  opt1.y_step = step_size;   opt1.x_win = 32; opt1.y_win = 32;
 opt2 = opt1;                       opt2.x_step = step_size;  opt2.y_step =step_size;    opt2.x_win = 32; opt2.y_win = 32;
 opt3 = opt1;                       opt3.x_step = step_size;  opt3.y_step = step_size;   opt3.x_win = 32; opt3.y_win = 32;
 pass = 3;
 tic;
 [x_WIDIM,y_WIDIM,u_WIDIM,v_WIDIM] =  PIV_analysis(I1,I2,pass,opt1,opt2,opt3);
 WIDIM_time = toc;
%- Parameters meaning can be found in <pivlab/Accuracy.m>
% s = {'Int. area 1',32;'Step size 1',16;'Subpix. finder',1;'Mask',[];'ROI',[];'Nr. of passes',3;'Int. area 2',32;'Int. area 3',16;'Int. area 4',32;'Window deformation','*spline'};
% p = {'ROI',s{5,2};'CLAHE',1;'CLAHE size',50;'Highpass',0;'Highpass size',15;'Clipping',0; 'Wiener',0;'Wiener size',3};
% image1 = PIVlab_preproc (I1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2}); %preprocess images
% image2 = PIVlab_preproc (I2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2});
% [y_WIDIM, x_WIDIM, v_WIDIM, u_WIDIM, ~] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});




%- Process the images with One-pass FFT-CC
 opt1.x_start = 5;opt1.y_start = 5; opt1.x_step = step_size; opt1.y_step = step_size;  opt1.x_win = 32; opt1.y_win = 32;
 opt2 = opt1; opt3 = opt1;                      
 pass = 1;
 tic
 [x_FFTCC,y_FFTCC,u_FFTCC,v_FFTCC] =  PIV_analysis(I1,I2,pass,opt1,opt2,opt3);
 FFTCC_time =toc;
%- Parameters meaning can be found in <pivlab/Accuracy.m>
% s = {'Int. area 1',32;'Step size 1',16;'Subpix. finder',1;'Mask',[];'ROI',[];'Nr. of passes',1;'Int. area 2',32;'Int. area 3',32;'Int. area 4',32;'Window deformation','*spline'};
% p = {'ROI',s{5,2};'CLAHE',1;'CLAHE size',50;'Highpass',0;'Highpass size',15;'Clipping',0; 'Wiener',0;'Wiener size',3};
% image1 = PIVlab_preproc (I1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2}); %preprocess images
% image2 = PIVlab_preproc (I2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2});
% [y_FFTCC, x_FFTCC, v_FFTCC, u_FFTCC, ~] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});


%- Process the images with Coarse-to-fine optical flow using Lucas&Kanade method.
[v_LKOF,u_LKOF] = opticalFlow(I1,I2,'smooth',1,'radius',16,'maxScale',1,'minScale',1/16,'type','LK');
[y_LKOF,x_LKOF] = meshgrid(1:size(v_LKOF,2),1:size(v_LKOF,1));


%- Process the images with pivdnn
opt.display = true; opt.stepSize = [step_size,step_size];opts.minPosi  = [4,4];opt.GPU = false;opt.Cudnn = false;
tic;
[x_pivdnn,y_pivdnn,u_pivdnn,v_pivdnn] =pivdnn(I1,I2,opt);
pivdnn_time = toc;
run_time = [FFTCC_time,WIDIM_time,pivdnn_time]

% u_pivdnn = smoothn(u_pivdnn,'robust');
% v_pivdnn = smoothn(v_pivdnn,'robust');

%- Resample the vector fields with the same nodes of pivdnn (step size = 16)
v_WIDIM = interp2(y_WIDIM,x_WIDIM,v_WIDIM,y_pivdnn,x_pivdnn,'nearest');
u_WIDIM = interp2(y_WIDIM,x_WIDIM,u_WIDIM,y_pivdnn,x_pivdnn,'nearest');
v_FFTCC = interp2(y_FFTCC,x_FFTCC,v_FFTCC,y_pivdnn,x_pivdnn,'nearest');
u_FFTCC = interp2(y_FFTCC,x_FFTCC,u_FFTCC,y_pivdnn,x_pivdnn,'nearest');
v_LKOF = interp2(y_LKOF,x_LKOF,v_LKOF,y_pivdnn,x_pivdnn,'nearest');
u_LKOF = interp2(y_LKOF,x_LKOF,u_LKOF,y_pivdnn,x_pivdnn,'nearest');

v_WIDIM(isnan(v_WIDIM)) = 0; u_WIDIM(isnan(u_WIDIM)) = 0; 
v_FFTCC(isnan(v_FFTCC)) = 0; u_FFTCC(isnan(u_FFTCC)) = 0; 
v_LKOF(isnan(v_LKOF)) = 0;   u_LKOF(isnan(u_LKOF)) = 0; 
v_pivdnn(isnan(v_pivdnn)) = 0;   u_pivdnn(isnan(u_pivdnn)) = 0; 

%- Display the results
if opts.show
    H0 = figure; % images
    set(gcf,'Name','The input PIV images');
    subplot 121; imshow(I1,[]);title('The first image');
    subplot 122; imshow(I2,[]);title('The second image');
    
    H1 = figure;% visualization of the vector field
    set(gcf,'Name','The PIV vectors');
    subplot 221;quiver(y_pivdnn(:),x_pivdnn(:),v_FFTCC(:),u_FFTCC(:),2);set(gca,'YDir','reverse'); title('One Pass FFT cross-correlation(FFTCC) PIV');
    subplot 222;quiver(y_pivdnn(:),x_pivdnn(:),v_WIDIM(:),u_WIDIM(:),2);set(gca,'YDir','reverse'); title('WIndow Deformation Iterative Multigrid (WIDIM) PIV');
    subplot 223;quiver(y_pivdnn(:),x_pivdnn(:),v_LKOF(:),u_LKOF(:),2);set(gca,'YDir','reverse'); title('Coarse-to-fine LK optical flow');
    subplot 224;quiver(y_pivdnn(:),x_pivdnn(:),v_pivdnn(:),u_pivdnn(:),2);set(gca,'YDir','reverse'); title('PIV-DCNN');
    
    H2 = figure;% contourf plot
    set(gcf,'Name','The PIV vector amplitude contours');
    subplot 221;set(gca,'YDir','reverse'); title('One Pass FFT cross-correlation(FFTCC) PIV');
    hold on; contourf(y_pivdnn,x_pivdnn,sqrt(v_FFTCC.^2+u_FFTCC.^2),[0.6 1 2 3 4 5 6 7 8],'b-');
    subplot 222;set(gca,'YDir','reverse'); title('WIndow Deformation Iterative Multigrid (WIDIM) PIV');
    hold on; contourf(y_pivdnn,x_pivdnn,sqrt(v_WIDIM.^2+u_WIDIM.^2),[0.6 1 2 3 4 5 6 7 8],'b-');
    subplot 223;set(gca,'YDir','reverse'); title('Coarse-to-fine LK optical flow');
    hold on; contourf(y_pivdnn,x_pivdnn,sqrt(v_LKOF.^2+u_LKOF.^2),[0.6 1 2 3 4 5 6 7 8],'b-');
    subplot 224;set(gca,'YDir','reverse'); title('PIV-DCNN');
    hold on; contourf(y_pivdnn,x_pivdnn,sqrt(v_pivdnn.^2+u_pivdnn.^2),[0.6 1 2 3 4 5 6 7 8],'b-');
    
    H3 = figure;% Vector component Histogram
    set(gcf,'Name','The PIV vector histogram');
    subplot 221;xlabel('\fontsize{14}Component displacement[pixel]');ylabel('\fontsize{14}Vector count');
    hist([v_FFTCC(:);u_FFTCC(:)],[-15:0.1:15]); title('One Pass FFT cross-correlation(FFTCC) PIV');
    subplot 222;xlabel('\fontsize{14}Component displacement[pixel]');ylabel('\fontsize{14}Vector count');
    hist([v_WIDIM(:);u_WIDIM(:)],[-15:0.1:15]); title('WIndow Deformation Iterative Multigrid (WIDIM) PIV');
    subplot 223;xlabel('\fontsize{14}Component displacement[pixel]');ylabel('\fontsize{14}Vector count');
    hist([v_LKOF(:);u_LKOF(:)],[-15:0.1:15]); title('Coarse-to-fine LK optical flow');
    subplot 224;xlabel('\fontsize{14}Component displacement[pixel]');ylabel('\fontsize{14}Vector count');
    hist([v_pivdnn(:);u_pivdnn(:)],[-15:0.1:15]); title('PIV-DCNN');
    
    sx = max(size(log(turb_energy_spectrum(u_FFTCC,v_FFTCC)+1)));
    w =1:sx;
    H4 = figure;% Turbulent Energy Spectrum
    set(gcf,'Name','The turbulent energy spectrum of vector field');
    grid off;hold on;box on;xlim([1-0.5 sx+0.5]); xlabel('\fontsize{14} Wave Number'); ylabel('\fontsize{14} (Energy+1) / Wave Number')
    title('\fontsize{18}Turbulent energy spectrum comparison')
    plot(w,(turb_energy_spectrum(u_FFTCC,v_FFTCC)+1),'-^','LineWidth',2,'Color',[0.0 0.0 1.0], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.0 0.0 1.0]);%original flow
    plot(w,(turb_energy_spectrum(u_WIDIM,v_WIDIM)+1),'--v','LineWidth',2,'Color',[0,0.45,0.74], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0,0.45,0.74]);%corrupted flow
    plot(w,(turb_energy_spectrum(u_LKOF,v_WIDIM)+1),'--h','LineWidth',2,'Color',[0.93,0.69,0.13], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.93,0.69,0.13]);%conventional methods without smoothness
    plot(w,(turb_energy_spectrum(u_pivdnn,v_pivdnn)+1),'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);%conventional methods
    H11 = legend('\fontsize{14}One Pass FFT cross-correlation(FFTCC) PIV','\fontsize{14}WIndow Deformation Iterative Multigrid (WIDIM) PIV','\fontsize{14}Coarse-to-fine LK optical flow','\fontsize{14}PIV-DCNN');
    set(gca,'Yscale','log'); set(gca,'fontsize',12)
    set(H4,'position',[100, 100, 800, 500]); 
    
    set(H1, 'position', get(0,'ScreenSize'));set(H2, 'position', get(0,'ScreenSize'));set(H3, 'position', get(0,'ScreenSize'));
    
    % output the data to origin, to generate better visualization

        Fig7a = [y_pivdnn(:),x_pivdnn(:),y_pivdnn(:)+8*v_FFTCC(:),x_pivdnn(:)+8*u_FFTCC(:),sqrt(v_FFTCC(:).^2+u_FFTCC(:).^2)];
        Fig7b = [y_pivdnn(:),x_pivdnn(:),y_pivdnn(:)+8*v_WIDIM(:),x_pivdnn(:)+8*u_WIDIM(:),sqrt(v_WIDIM(:).^2+u_WIDIM(:).^2)];
        Fig7c = [y_pivdnn(:),x_pivdnn(:),y_pivdnn(:)+8*v_pivdnn(:),x_pivdnn(:)+8*u_pivdnn(:),sqrt(v_pivdnn(:).^2+u_pivdnn(:).^2)];
    
%     v_pivdnn(v_pivdnn>3) = 3; v_pivdnn(v_pivdnn<-3) = -3;
%     u_pivdnn(u_pivdnn>3) = 3; u_pivdnn(u_pivdnn<-3) = -3;
%     v_WIDIM(v_WIDIM>3) = 3; v_WIDIM(v_WIDIM<-3) = -3;
%     u_WIDIM(u_WIDIM>3) = 3; u_WIDIM(u_WIDIM<-3) = -3;
%     v_FFTCC(v_FFTCC>3) = 3; v_FFTCC(v_FFTCC<-3) = -3;
%     u_FFTCC(u_FFTCC>3) = 3; u_FFTCC(u_FFTCC<-3) = -3;
%     %- Fig11a
%     figure; imagesc(-u_pivdnn); set(gca,'ydir','reverse');colormap('jet');set(gca,'ydir','reverse');colormap('jet'); set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
%     %- Fig11b
%     figure; imagesc(-[u_FFTCC(115:150,5:32),-ones(36,1),u_WIDIM(115:150,5:32),-ones(36,1),u_pivdnn(115:150,5:32);...
%                       -ones(1,86);
%                       u_FFTCC(160:190,130:157),-ones(31,1),u_WIDIM(160:190,130:157),-ones(31,1),u_pivdnn(160:190,130:157)]); set(gca,'ydir','reverse');colormap('jet'); set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
% 

end

end

