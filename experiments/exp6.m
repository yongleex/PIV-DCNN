%% Experiments 6
%- Test the frequence response of PIV-DNN (modulation)
%- Yong Lee (leeyong@hust.edu.cn)
%- 2017-05-27

function exp6()
close all; clear

%- Add methods pathes
addpath('../PIV_DNN/','../PIVCC/','../pivlab/','../OpticalFlow/');

%% Initial settings
NormalFrequence = 0.0:0.1:1.6; %NF = Win_size/wave_length;

%% Main Loop for different scale
for j = 1:10
    for i = 1:numel(NormalFrequence)
        [j,i]
        %% moving average filter with size 32
        %- Synthesize the image pair
        Win_size = 32;
        wave_length = Win_size/NormalFrequence(i);
        
        opts.size = [128,128]; % Image size
        opts.u =  @(x,y) 1.5*sin(2*pi*y/wave_length+pi/2);
        [xg,yg,u,v,I1,I2] = synSimu(opts);
        %     imwrite(uint8(I1),'I1.tif');
        %     imwrite(uint8(I2),'I2.tif');
        
        %- Moving average filter
        weights = ones(size(u));
        Filter  = ones(Win_size);
        u_MA = conv2(u,Filter,'same')./conv2(weights,Filter,'same');
        M_MA(i,j) = modul(u,u_MA);
        
        
        %- Process the images with One-pass FFT-CC
        step_size= 16;
        opt1.x_start = 5;opt1.y_start = 5; opt1.x_step = step_size; opt1.y_step = step_size;  opt1.x_win = 32; opt1.y_win = 32;
        opt2 = opt1; opt3 = opt1;    pass = 1;
        [x_FFTCC,y_FFTCC,u_FFTCC,v_FFTCC] =  PIV_analysis(I1,I2,pass,opt1,opt2,opt3);
        v_truth = interp2(yg,xg,v,y_FFTCC,x_FFTCC,'nearest');
        u_truth = interp2(yg,xg,u,y_FFTCC,x_FFTCC,'nearest');
        M_FFTCC(i,j) = modul(u_truth,u_FFTCC);
        
        Win_size = 32; step_size = 16;
        wave_length = Win_size/NormalFrequence(i);
        opts.u =  @(x,y) 1.5*sin(2*pi*y/wave_length+pi/2);
        [xg,yg,u,v,I1,I2] = synSimu(opts);
        
        opt1.x_start = 5;opt1.y_start = 5; opt1.x_step = 16; opt1.y_step = 16;   opt1.x_win = 32; opt1.y_win = 32;
        opt2 = opt1;                       opt2.x_step = 16; opt2.y_step = 16;   opt2.x_win = 32; opt2.y_win = 32;
        opt3 = opt1;                       opt3.x_step = 16; opt3.y_step = 16;   opt3.x_win = 32; opt3.y_win = 32;
        pass = 3;
        [x_WIDIM,y_WIDIM,u_WIDIM,v_WIDIM] =  PIV_analysis(I1,I2,pass,opt1,opt2,opt3);
        
%         %- Parameters meaning can be found in <pivlab/Accuracy.m>
%         s = {'Int. area 1',32;'Step size 1',16;'Subpix. finder',1;'Mask',[];'ROI',[];'Nr. of passes',1;'Int. area 2',32;'Int. area 3',16;'Int. area 4',16;'Window deformation','*spline'};
%         p = {'ROI',s{5,2};'CLAHE',1;'CLAHE size',50;'Highpass',0;'Highpass size',15;'Clipping',0; 'Wiener',0;'Wiener size',3};
%         image1 = PIVlab_preproc (I1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2}); %preprocess images
%         image2 = PIVlab_preproc (I2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2});
%         [y_WIDIM, x_WIDIM, v_WIDIM, u_WIDIM, ~] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
        
        v_truth = interp2(yg,xg,v,y_WIDIM,x_WIDIM,'nearest');
        u_truth = interp2(yg,xg,u,y_WIDIM,x_WIDIM,'nearest');
%         u_WIDIM = smoothn(u_WIDIM,'robust');
        M_WIDIM(i,j) = modul(u_truth,u_WIDIM);
        
        %- Process the images with pivdnn
        Win_size = 64; step_size = 16;
        wave_length = Win_size/NormalFrequence(i);
        opts.u =  @(x,y) 1.5*sin(2*pi*y/wave_length+pi/2);
        [xg,yg,u,v,I1,I2] = synSimu(opts);
        
        opt.display = false; opt.stepSize = [step_size,step_size];opts.minPosi  = [4,4];opt.GPU = false;opt.Cudnn = false;
        [x_pivdnn,y_pivdnn,u_pivdnn,v_pivdnn] =pivdnn(I1,I2,opt);
        
        v_truth = interp2(yg,xg,v,y_pivdnn,x_pivdnn,'nearest');
        u_truth = interp2(yg,xg,u,y_pivdnn,x_pivdnn,'nearest');
        M_pivdnn(i,j) = modul(u_truth,u_pivdnn);
    end
end

save exp6
M_MA = nanmean(M_MA,2);
M_FFTCC = nanmean(M_FFTCC,2);
M_WIDIM = nanmean(M_WIDIM,2);
M_pivdnn = nanmean(M_pivdnn,2);

%% Display the results
Colours = linspecer(5);
lineStyle = {'-'; '-.'; 'd';'p'; '--o';};
H = figure; hold on; grid on; box on;
xlim([min(NormalFrequence(:)),max(NormalFrequence(:))]);
ylim([-0.5,1.2]);
set(gca,'ytick',[-0.5 0 0.5 0.8 1 1.2]);
plot(NormalFrequence,sinc(NormalFrequence),lineStyle{1},'LineWidth',2,'color',Colours(1,:));
plot(NormalFrequence,M_MA,lineStyle{2},'LineWidth',2,'color',Colours(2,:))
plot(NormalFrequence,M_FFTCC,lineStyle{3},'LineWidth',2,'MarkerFaceColor',Colours(3,:), 'MarkerEdgeColor','k')
plot(NormalFrequence,M_WIDIM,lineStyle{4},'LineWidth',2,'MarkerFaceColor',Colours(4,:), 'MarkerEdgeColor','k')
plot(NormalFrequence,M_pivdnn,lineStyle{5},'LineWidth',2,'color',Colours(5,:),'MarkerFaceColor',Colours(5,:), 'MarkerEdgeColor','k')
H11 = legend('\fontsize{16}sinc(theory)','\fontsize{16}MA Filter','\fontsize{16}FFTCC','\fontsize{16}WIDIM','\fontsize{16}PIV-DCNN','Location','SouthWest');
xlabel('\fontsize{16}Normalized Window Size');
ylabel('\fontsize{16}Modulation Coefficient');set(gca,'fontsize',16); set(H,'position',[400 200 600 350]);
plot(NormalFrequence,0.9*ones(size(NormalFrequence)),'--b','LineWidth',2);


end

%% Modulation coefficient calculation
function  m = modul(u,u_cal)
m = (sum(u_cal(:).*u(:))/(sum(u(:).^2))); % Eq. x
end


%% Generate random artificial particle images
function [xg,yg,u,v,I1,I2] = synSimu(opts)
%- Default settings
opt.size = [512,512]; % Image size
opt.Intens  = 0.8;  % particle intensity
opt.dIntens = 0.001; % particle intensity variation


opt.ppp  = 0.025;%particle concentration (particles per pixel)
opt.dt   = 2.2; %particle diameter
opt.ddt  = 0.3;   %particle diameter variation
opt.inPlaneRatio = 0.9999; % The ratio of the in-plane particles

ud = @(x,y)  2 + 0*x;  % default uniform displacement
vd = @(x,y)  0 + 0*x;
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


