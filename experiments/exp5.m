%% Experiments 5 (comparison with other PIV evaluation methods, uniform flow)
%- Compare the PIV-DNN performance with state-of-art PIV methods in the synthetic uniform flow cases.
%- Adopt the root mean square (RMS) error
%- Yong Lee (leeyong@hust.edu.cn)
%- 2016-12-18
%- Modified @ 2017-06-02
function exp5()
close all; clear;
%- Configuration ()
opts.MonteCarlo = 100; % 100 to run the figure in our manuscript
opts.show = false; % you can set it as true with debug model to see the intermediate results

opts.size = [128,128]; % Image size
opts.ppp  = 0.05;%particle concentration (particles per pixel)
opts.dt   = 2.2; %particle diameter
opts.inPlaneRatio = 0.999; % The ratio of the in-plane particles

%- Add methods pathes
addpath('../PIV_DNN/','../PIVCC/','../pivlab/','../OpticalFlow/');

Displacement = 0:0.1:10;
BiasErr = zeros(3,numel(Displacement)); RMSErr = zeros(3,numel(Displacement));
for j = 1:numel(Displacement)
    pivdnn_Res = []; WIDIM_Res = [];  FFTCC_Res = []; LKOF_Res = [];
    for i = 1:opts.MonteCarlo
        [j,i]
        %- Generating the synthetic vector field and images
        ud = @(x,y)  Displacement(j) + 0*x;
        vd = @(x,y)  0 + 0*x;
        opts.u = ud; opts.v = vd;
        [x,y,u,v,I1,I2] = synSimu(opts);
        
        %- Add noise to the images
        %         I1 = imnoise(I1,'gaussian',0,0.001);
        %         I2 = imnoise(I2,'gaussian',0,0.001);
        
        %- Process the images with One-pass FFT-CC
        opt1.x_start = 5;opt1.y_start = 5; opt1.x_step = 16; opt1.y_step = 16;  opt1.x_win = 32; opt1.y_win = 32;
        opt2 = opt1; opt3 = opt1;
        pass = 1;
        [x_FFTCC,y_FFTCC,u_FFTCC,v_FFTCC] =  PIV_analysis(I1,I2,pass,opt1,opt2,opt3);
        
        %- Parameters meaning can be found in <pivlab/Accuracy.m>
        % s = {'Int. area 1',32;'Step size 1',16;'Subpix. finder',1;'Mask',[];'ROI',[];'Nr. of passes',1;'Int. area 2',32;'Int. area 3',32;'Int. area 4',32;'Window deformation','*spline'};
        % p = {'ROI',s{5,2};'CLAHE',1;'CLAHE size',50;'Highpass',0;'Highpass size',15;'Clipping',0; 'Wiener',0;'Wiener size',3};
        % image1 = PIVlab_preproc (I1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2}); %preprocess images
        % image2 = PIVlab_preproc (I2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2});
        % [y_FFTCC, x_FFTCC, v_FFTCC, u_FFTCC, ~] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
        
        %- Process the images with WIndow Deformation Iterative Multigrid (WIDIM) PIV
        opt1.x_start = 5;opt1.y_start = 5; opt1.x_step = 16; opt1.y_step = 16;  opt1.x_win = 32; opt1.y_win = 32;
        opt2 = opt1;                       opt2.x_step = 8;  opt2.y_step = 8;   opt2.x_win = 16; opt2.y_win = 16;
        opt3 = opt1;                       opt3.x_step = 16; opt3.y_step = 16;   opt3.x_win = 16; opt3.y_win = 16;
        pass = 3;
        [x_WIDIM,y_WIDIM,u_WIDIM,v_WIDIM] =  PIV_analysis(I1,I2,pass,opt1,opt2,opt3);
        
        %- Parameters meaning can be found in <pivlab/Accuracy.m>
        % s = {'Int. area 1',32;'Step size 1',16;'Subpix. finder',1;'Mask',[];'ROI',[];'Nr. of passes',3;'Int. area 2',32;'Int. area 3',16;'Int. area 4',32;'Window deformation','*spline'};
        % p = {'ROI',s{5,2};'CLAHE',1;'CLAHE size',50;'Highpass',0;'Highpass size',15;'Clipping',0; 'Wiener',0;'Wiener size',3};
        % image1 = PIVlab_preproc (I1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2}); %preprocess images
        % image2 = PIVlab_preproc (I2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2});
        % [y_WIDIM, x_WIDIM, v_WIDIM, u_WIDIM, ~] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
        
        %- Process the images with Coarse-to-fine optical flow using Lucas&Kanade method.
        % [v_HSOF,u_HSOF] = opticalFlow(I1,I2,'smooth',1,'radius',10,'maxScale',1,'minScale',1/16,'type','HS');
        [v_LKOF,u_LKOF] = opticalFlow(I1,I2,'smooth',1,'radius',16,'maxScale',1,'minScale',1/16,'type','LK');
        [y_LKOF,x_LKOF] = meshgrid(1:size(v_LKOF,2),1:size(v_LKOF,1));
        
        
        %- Process the images with pivdnn
        opt.display = false; opt.stepSize = [16,16];opts.minPosi  = [4,4];opt.GPU = false;opt.Cudnn = false;
        [x_pivdnn,y_pivdnn,u_pivdnn,v_pivdnn] =pivdnn(I1,I2,opt);
        
        %- Resample the vector fields with the same nodes of pivdnn (step size = 16)
        u = interp2(y,x,u,y_pivdnn,x_pivdnn,'nearest');
        v = interp2(y,x,v,y_pivdnn,x_pivdnn,'nearest');
        v_WIDIM = interp2(y_WIDIM,x_WIDIM,v_WIDIM,y_pivdnn,x_pivdnn,'nearest');
        u_WIDIM = interp2(y_WIDIM,x_WIDIM,u_WIDIM,y_pivdnn,x_pivdnn,'nearest');
        v_FFTCC = interp2(y_FFTCC,x_FFTCC,v_FFTCC,y_pivdnn,x_pivdnn,'nearest');
        u_FFTCC = interp2(y_FFTCC,x_FFTCC,u_FFTCC,y_pivdnn,x_pivdnn,'nearest');
        v_LKOF = interp2(y_LKOF,x_LKOF,v_LKOF,y_pivdnn,x_pivdnn,'nearest');
        u_LKOF = interp2(y_LKOF,x_LKOF,u_LKOF,y_pivdnn,x_pivdnn,'nearest');
        
        %- Display the intermediate results (images, velocity fields ground truth and estimations)
        if opts.show
            figure;
            subplot 221; imshow(I1,[]); title('Image 1');
            subplot 222; imshow(I2,[]); title('Image 2');
            subplot 223;quiver(y_pivdnn(:),x_pivdnn(:),v(:),u(:));set(gca,'YDir','reverse'); title('Truth');
            figure;
            subplot 221;quiver(y_pivdnn(:),x_pivdnn(:),v_FFTCC(:),u_FFTCC(:),2);set(gca,'YDir','reverse'); title('One Pass FFT cross-correlation(FFTCC) PIV');
            subplot 222;quiver(y_pivdnn(:),x_pivdnn(:),v_WIDIM(:),u_WIDIM(:),2);set(gca,'YDir','reverse'); title('WIndow Deformation Iterative Multigrid (WIDIM) PIV');
            subplot 223;quiver(y_pivdnn(:),x_pivdnn(:),v_LKOF(:),u_LKOF(:),2);set(gca,'YDir','reverse'); title('Coarse-to-fine HS optical flow');
            subplot 224;quiver(y_pivdnn(:),x_pivdnn(:),v_pivdnn(:),u_pivdnn(:),2);set(gca,'YDir','reverse'); title('pivdnn');
        end
        
        pivdnn_Res  = [pivdnn_Res; (u_pivdnn(:)-Displacement(j))];
        WIDIM_Res   = [WIDIM_Res;  (u_WIDIM(:) -Displacement(j))];
        FFTCC_Res   = [FFTCC_Res;  (u_FFTCC(:) -Displacement(j))];
        LKOF_Res    = [LKOF_Res;   (u_LKOF(:)  -Displacement(j))];
    end
    %- Set the error to 1 when confront the outliers
    WIDIM_Res(isnan(WIDIM_Res))= 1;     FFTCC_Res(isnan(FFTCC_Res))= 1;
    
    
    %- Analysis the errors at a specified velocity.
    BiasErr(1,j) = mean(LKOF_Res); BiasErr(2,j) = mean(FFTCC_Res); BiasErr(3,j) = mean(WIDIM_Res); BiasErr(4,j) = mean(pivdnn_Res);
    RMSErr(1,j)  = sqrt(mean(LKOF_Res(:).^2)); RMSErr(2,j)  = sqrt(mean(FFTCC_Res(:).^2)); RMSErr(3,j)  = sqrt(mean(WIDIM_Res(:).^2)); RMSErr(4,j)  = sqrt(mean(pivdnn_Res(:).^2));
end
save exp5; % save the results

%- plot the results(bias error and RMS error)
H1 = figure; hold on; box on; %title('bias error');
xlim([min(Displacement), max(Displacement)]);%ylim([-0.01 1.01]);set(gca,'yscale','log');
plot(Displacement,BiasErr(1,:),'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);
plot(Displacement,BiasErr(2,:),'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
plot(Displacement,BiasErr(3,:),'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(Displacement,BiasErr(4,:),'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
H11 = legend('\fontsize{16}LKOF','\fontsize{16}FFTCC','\fontsize{16}WIDIM','\fontsize{16}PIV-DCNN');
xlabel('\fontsize{16}Particle image displacement \deltax[pixel]');
ylabel('\fontsize{16}Mean error(Bias) [pixel] ');set(gca,'fontsize',16); set(H1,'position',[300 100 600 350]);

H2 = figure; hold on; box on; %title('rms error');
xlim([min(Displacement), max(Displacement)]);%ylim([-0.01 1.01]);set(gca,'yscale','log');
plot(Displacement,RMSErr(1,:),'--+','LineWidth',2,'Color',[0.3,0.75,0.93], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3,0.75,0.93]);
plot(Displacement,RMSErr(2,:),'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
plot(Displacement,RMSErr(3,:),'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(Displacement,RMSErr(4,:),'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
H11 = legend('\fontsize{16}LKOF','\fontsize{16}FFTCC','\fontsize{16}WIDIM','\fontsize{16}PIV-DCNN');
xlabel('\fontsize{16}Particle image displacement \deltax[pixel]');
ylabel('\fontsize{16}RMS error [pixel] ');set(gca,'fontsize',16); set(H2,'position',[300 100 600 350]);
%- Bias without LKOF
H3 = figure; hold on; box on; %title('bias error');
xlim([min(Displacement), max(Displacement)]);%ylim([-0.01 1.01]);set(gca,'yscale','log');
plot(Displacement,BiasErr(2,:),'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
plot(Displacement,BiasErr(3,:),'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(Displacement,BiasErr(4,:),'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
H11 = legend('\fontsize{16}FFTCC','\fontsize{16}WIDIM','\fontsize{16}PIV-DCNN');
xlabel('\fontsize{16}Particle image displacement \deltax[pixel]');
ylabel('\fontsize{16}Mean error(Bias) [pixel] ');set(gca,'fontsize',16); set(H3,'position',[300 100 600 350]);
%- RMS without LKOF
H4 = figure; hold on; box on; grid on;%title('rms error');
xlim([min(Displacement), max(Displacement)]);ylim([0.0001 0.08]); %set(gca,'yscale','log');

plot(Displacement,RMSErr(2,:),'--d','LineWidth',2,'Color',[0.85,0.33,0.1], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85,0.33,0.1]);
plot(Displacement,RMSErr(3,:),'--p','LineWidth',2,'Color',[0.64,0.08,0.18], 'MarkerSize',6, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.64,0.08,0.18]);
plot(Displacement,RMSErr(4,:),'-o','LineWidth',2,'Color',[0.49,0.18,0.56], 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.49,0.18,0.56]);
H11 = legend('\fontsize{16}FFTCC','\fontsize{16}WIDIM','\fontsize{16}PIV-DCNN');
xlabel('\fontsize{16}Particle image displacement \deltax[pixel]');
ylabel('\fontsize{16}RMS error [pixel] ');set(gca,'fontsize',16); set(H4,'position',[300 100 600 350]);

end

%% Generate random artificial particle images
function [xg,yg,u,v,I1,I2] = synSimu(opts)
%- Default settings
opt.size = [128,128]; % Image size
opt.Intens  = 0.8;  % particle intensity
opt.dIntens = 0.001; % particle intensity variation


opt.ppp  = 0.01;%particle concentration (particles per pixel)
opt.dt   = 4; %particle diameter
opt.ddt  = 0.3;   %particle diameter variation
opt.inPlaneRatio = 0.9999; % The ratio of the in-plane particles

ud = @(x,y)  2 + 0*x;  % default uniform displacement
vd = @(x,y) -1 + 0*x;
opt.u = ud; opt.v = vd;

%- User settings
opts.name = [];
opt = argparse(opt, opts); % using the usr settings

%- Initial calculation
partAm = round(opt.ppp*prod(opt.size)); % particle number
[yg,xg] = meshgrid(1:opt.size(1),1:opt.size(2));
u = opt.u(xg,yg); v = opt.v(xg,yg);

%- Generate the particles
Particle.position(1,:) = opt.size(1)*rand(1,partAm,'single');
Particle.position(2,:) = opt.size(2)*rand(1,partAm,'single');
Particle.diameter  = opt.dt + opt.ddt* randn(1,partAm,'single');
Particle.intensity = opt.Intens+ opt.dIntens * randn(1,partAm,'single');
Particle.number = floor(partAm*opt.inPlaneRatio);

Particle1 = Particle;
Particle2 = Particle; % movement
Particle2.position(1,:) = Particle2.position(1,:) + opt.u(Particle2.position(1,:),Particle2.position(2,:));
Particle2.position(2,:) = Particle2.position(2,:) + opt.v(Particle2.position(1,:),Particle2.position(2,:));

%- remove the out-of-Plane particles
Index1 = randperm(partAm); Index1(1:round(partAm*opt.inPlaneRatio)) = [];
Index2 = randperm(partAm); Index2(1:round(partAm*opt.inPlaneRatio)) = [];
Particle1.position(:,Index1) = [];Particle1.diameter(Index1) = [];Particle1.intensity(Index1) = [];
Particle2.position(:,Index2) = [];Particle2.diameter(Index2) = [];Particle2.intensity(Index2) = [];

%- Generate the images
I1 = genImage(Particle1,xg,yg,opt); I1 = uint8(I1*255);
I2 = genImage(Particle2,xg,yg,opt); I2 = uint8(I2*255);
end

%- Copy the field in opt to opts
function opts = argparse(opts, opt)
opts = opts;
fieldName = fieldnames(opt);
for i = 1:numel(fieldName)
    eval(cell2mat(['opts.',fieldName(i),' =', 'opt.',fieldName(i),';']))
end
end

%% Generate the images from particles
function I1 = genImage(Particle,xg,yg,opt)
I1 = zeros(opt.size,'single');

% % imaging procedure for each particle (Gaussian theory but very slow)
% for i=1:1:Particle.number
%     I1  = I1+ (Particle.intensity(i).*exp(-8*((xg-Particle.position(1,i)).^2+(yg-Particle.position(2,i)).^2)./(Particle.diameter(i)^2)));
% end

%- imaging procedure for each particle
Idx1 = round(max( Particle.position(1,:)-Particle.diameter,1));
Idx2 = round(min( Particle.position(1,:)+Particle.diameter,opt.size(1)));
Idy1 = round(max( Particle.position(2,:)-Particle.diameter,1));
Idy2 = round(min( Particle.position(2,:)+Particle.diameter,opt.size(2)));

for i =1:Particle.number
    %     I1(Idx1(i):Idx2(i),Idy1(i):Idy2(i)) = I1(Idx1(i):Idx2(i),Idy1(i):Idy2(i)) + ...
    %         (Particle.intensity(i).*exp(-8*((xg(Idx1(i):Idx2(i),Idy1(i):Idy2(i))-Particle.position(1,i)).^2+(yg(Idx1(i):Idx2(i),Idy1(i):Idy2(i))-Particle.position(2,i)).^2)./(Particle.diameter(i)^2)));
    
    dr = Particle.diameter(i);
    I1(Idx1(i):Idx2(i),Idy1(i):Idy2(i)) = I1(Idx1(i):Idx2(i),Idy1(i):Idy2(i)) + ...
        Particle.intensity(i)./(erf(2.8284*(0.5)/dr)).* ...
        (erf(2.8284*((xg(Idx1(i):Idx2(i),Idy1(i):Idy2(i))-Particle.position(1,i))+0.5)/dr)-erf(2.8284*((xg(Idx1(i):Idx2(i),Idy1(i):Idy2(i))-Particle.position(1,i))-0.5)/dr)).*(erf(2.8284*((yg(Idx1(i):Idx2(i),Idy1(i):Idy2(i))-Particle.position(2,i))+0.5)/dr)-erf(2.8284*((yg(Idx1(i):Idx2(i),Idy1(i):Idy2(i))-Particle.position(2,i))-0.5)/dr))/4;
end

I1 = single(I1);
end