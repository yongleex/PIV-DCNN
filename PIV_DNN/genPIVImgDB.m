%% Generate the synthetic PIV images with given number
% (down) x-direction with first component u
% (right) y-direction with second vector component v
function [images, meta] = genPIVImgDB(DataSize,opt)
% clear;
%- Check the input data size
if nargin >= 1,DatabaseSize = DataSize;else close all; DatabaseSize = 1;end

%- Init the size of images
opts.Ni = [95,95];    % Full image size with boundary motion problem
opts.No = [64,64];    % Output Image size by cropping

%- Given a range to the particle's features
opts.pdr = [0.35,6]; % particle diameter range in pixel
opts.pir = [0.8,1]; % particle intensity range
opts.pcr = [0.008,0.1]; % particle concentration/density (particles per pixel)
opts.por = [0,0.05];    % range of the percentage of the out-of-plane particles

%- Some configuration for vector field
opts.mvr = [0,5.0];     % maximal velocity component range.

opt.name = [];
opts = argparse(opts, opt); % using the usr settings

images.data   =  zeros(opts.No(1),opts.No(2),2,DatabaseSize,'single');% The output image patches
images.vector =  zeros(2,DatabaseSize,'single');        % The output displacement vectors with the same size of image patches

%- some pre-calculated value for a faster implementation
%- The grid nodes
[yg,xg]=meshgrid(linspace(1,opts.Ni(1),opts.Ni(1)),linspace(1,opts.Ni(2),opts.Ni(2)));

rng('shuffle'); % seeds the random number generator based on the current time
hwait=waitbar(0,'Generating start. Please wait>>>>>>>>');
for i = 1:DatabaseSize
    %- Progress bar/ monitor
    i
    if rem(i,100)== 0, waitbar(i/DatabaseSize,hwait,'Generating PIV image patches >>>>>>');end
    
    %- The velocity distribution
    P= 0.1*(2*rand(1,23)-1);
    P([2,3,4,9,10,11,15,16,17,18,19,20]) = 100*P([2,3,4,9,10,11,15,16,17,18,19,20]);
    P(1) = single(unifrnd(-opts.mvr(2),opts.mvr(2)));
    P(8) = single(unifrnd(-opts.mvr(2),opts.mvr(2)));
    P(5) = 32* P(5);     P(12) = 32* P(12); % The range below 3.2 pixel
    P(6) = 20*pi*P(6);   P(13) = 20*pi*P(13); % the wave lambda >64
    P(7) = 20*pi*P(7);   P(14) = 20*pi*P(14); % The wave phase shift(0,2pi)
    P(21) = 0.25* P(21); P(22) = 20*pi*P(22); P(23) = 20*pi*P(23);
    if rand < 0.05 % Four special simple case
        Tp = P; P = 0*P;
        switch randi(4),
            case 1, P(1) = Tp(1); % uniform flow along a direction
            case 2, P([1,8]) = Tp([1,8]);% uniform flow has two components
            case 3, P([1, 2]) = Tp([1, 2]); % linear
            case 4, P([1, 5, 6, 7]) = Tp([1, 5, 6, 7] ); % sine wave and
        end
    end
    P = single(P);
    u = @(x,y) P(1)+ P(2)*y+ P(3)*y.^2 + P(4)*y.^3  +  P(5)*sin(P(6)*y  +P(7))-P(5)*sin(P(7)) + P(15)*x + 0.5*P(16)*x.^2 +     P(17)*x.*y + ...
        P(18)*x.*x.*x/3 + P(19)*x.*y.*y + P(20)*x.*x.*y/2 + P(21).*sin(P(22)*x+P(23)).*cos(P(22)*y+P(23)) - P(21).*sin(P(23)).*cos(P(23));
    v = @(x,y) P(8)+ P(9)*x+P(10)*x.^2 + P(11)*x.^3 + P(12)*sin(P(13)*x+P(14))-P(12)*sin(P(14)) - P(15)*y -     P(16)*x.*y - 0.5*P(17)*y.^2 - ...
        P(18)*x.*x.*y - P(19)*y.^3/3  - P(20)*x.*y.*y/2 - P(21).*cos(P(22)*x+P(23)).*sin(P(22)*y+P(23))+ P(21).*cos(P(23)).*sin(P(23));
    %     u(0,0), v(0,0),
    if rand < 0.6 % synthetic the images from particles
        ppp = single(unifrnd(opts.pcr(1),opts.pcr(2)));    % a specified case with a random particle density: particles per pixel
        pip = single(1- unifrnd(opts.por(1),opts.por(2))); % Take the out-of-plane motion into account: particles in-plane-motion percentage
        
        %%%%%%%%%%%%%using the * denote the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NpT= 10 means **********  for example, take 10 particles%%%%%%%%%%%%%%%%%
        % Np = 8  means ********    for the first image %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 ********  for the seconde image %%%%%%%%%%%%%%%%%%%%%%%%%
        Np  = round(prod(opts.Ni)*ppp); % particle numbers in an image
        NpT = round(prod(opts.Ni)*ppp/(pip)^2);
        
        maxV = single(unifrnd(opts.mvr(1),opts.mvr(2)));
        
        
        %- All the particles information before deformation
        Particle.number = NpT;
        Particle.position  = single([unifrnd(1,opts.Ni(1),[1,NpT]);unifrnd(1,opts.Ni(2),[1,NpT])]);
        Particle.intensity = single(unifrnd(opts.pir(1),opts.pir(2),[1,NpT]));
        Particle.diameter  = single((unifrnd(opts.pdr(1),opts.pdr(2)))* ones([1,NpT]));
        Particle1 = Particle; Particle2 = Particle;
        %- Get the displacement for every particle 
        [Vec,Vx,Vy] = genVecField(Particle2.position,maxV,opts,u,v);
        Particle1.position = Particle1.position - 0.5*Vec;
        Particle2.position = Particle2.position + 0.5*Vec;
        
        %     single(unifrnd(opts.pdr(1),unifrnd(mean(opts.pdr),opts.pdr(2)),[1,NpT]));
        
        %- Get the imaging particles.
         Particle1.number = Np; Particle1.position(:,(Np+1):end) = []; Particle1.intensity((Np+1):end) = []; Particle1.diameter((Np+1):end) = [];
         Particle2.number = Np; Particle2.position(:,1:(end-Np)) = []; Particle2.intensity(1:(end-Np)) = []; Particle2.diameter(1:(end-Np)) = [];
        
        %- Generate the images from particles without noise
        I1 = genImage(Particle1,xg,yg,opts);
        I2 = genImage(Particle2,xg,yg,opts);
        
        %- Add some noise to the images
        I1 = AddNoise(I1,xg,yg,opts); I2 = AddNoise(I2,xg,yg,opts);
    else % get the first image from images
        Vx = u(0,0);
        Vy = v(0,0);
        [I1,I2] = clipAndDefor(u,v,xg,yg,opts);
        %- Add some noise to the images
        I1 = AddNoise(I1,xg,yg,opts); I2 = AddNoise(I2,xg,yg,opts);
    end

    
    
    
    %     imwrite(uint8(I1*255),['I',num2str(i),'a.tif']);
    %     imwrite(uint8(I2*255),['I',num2str(i),'b.tif']);
    
    %- image data and labels
    I1 = I1(((opts.Ni-opts.No+1)/2:(opts.Ni+opts.No-1)/2)+0,((opts.Ni-opts.No+1)/2:(opts.Ni+opts.No-1)/2)+0);
    images.data(:,:,1,i) = I1-mean(I1(:)); % substracting the mean value
    I2 = I2(((opts.Ni-opts.No+1)/2:(opts.Ni+opts.No-1)/2)+0,((opts.Ni-opts.No+1)/2:(opts.Ni+opts.No-1)/2)+0);
    images.data(:,:,2,i) = I2-mean(I2(:));
    images.vector(1:2,i) = [Vx;Vy];
    
    if DatabaseSize ==1 % For debug display & output for figure 3
        figure(); imshow(I1,[]); title('Image 1');
        figure(); imshow(I2,[]); title('Image 2');
        imwrite(uint8(I1*255),'I1.tif');
        imwrite(uint8(I2*255),'I2.tif');
        [y,x] = meshgrid(-0.3368:1/9:0.3368); Vx =u(x,y);Vy =v(x,y);std(Vx(:));std(Vy(:)); figure; quiver(y,x,Vy,Vx);title('Flow pattern');  set(gca, 'YDir','reverse'); drawnow;
        figure(); mesh(x,y,double(Vx));
        figure(); mesh(x,y,double(Vy));
    end
    
end
close(hwait);

%- image set category
images.set(1:DatabaseSize) =1;
images.set(round(0.833333*DatabaseSize):round(0.916667*DatabaseSize)) =2;
images.set(round(0.916667*DatabaseSize):DatabaseSize) =3;

%- image data mean
images.data_mean = zeros([64,64,2]);

meta.sets{1} = 'train';meta.sets{2} = 'val';meta.sets{3} = 'test';
meta.classes ={'1','2','3'};
end



%% Generate the images from particles
function I1 = genImage(Particle,xg,yg,opt)
I1 = zeros(opt.Ni,'single');

% Add some noise to the particles
% Particle.position  = Particle.position  + unifrnd(0,0.50)*rand(size(Particle.position));
Particle.intensity = Particle.intensity + unifrnd(0,0.03,size(Particle.intensity));
Particle.diameter  = Particle.diameter  + unifrnd(0,0.25,size(Particle.diameter));

% % imaging procedure for each particle (Gaussian theory but very slow)
% for i=1:1:Particle.number
%     I1  = I1+ (Particle.intensity(i).*exp(-8*((xg-Particle.position(1,i)).^2+(yg-Particle.position(2,i)).^2)./(Particle.diameter(i)^2)));
% end

%- imaging procedure for each particle
Idx1 = round(max( Particle.position(1,:)-Particle.diameter,1));
Idx2 = round(min( Particle.position(1,:)+Particle.diameter,opt.Ni(1)));
Idy1 = round(max( Particle.position(2,:)-Particle.diameter,1));
Idy2 = round(min( Particle.position(2,:)+Particle.diameter,opt.Ni(2)));

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


%% Get a random vector field by randomly sampling
function [VecField,Vx,Vy] = genVecField(Pos,maxV,opt,u,v)
VecField = zeros(size(Pos),'single');
Pos(1,:) = Pos(1,:)./opt.Ni(1)-0.5; Pos(2,:) = Pos(2,:)./opt.Ni(2)-0.5; % Normalized the positions into [0,1]

VecField(1,:) = u(Pos(1,:),Pos(2,:));
VecField(2,:) = v(Pos(1,:),Pos(2,:));
% quiver(Pos(1,:),Pos(2,:),VecField(1,:),VecField(2,:))

Vx = u(0,0);
Vy = v(0,0);
% coef = maxV./max(abs([Vx,Vy])); % Take the maximal constraint into account
% VecField = coef*VecField; Vx = coef*Vx; Vy = coef*Vy;

% %- First flow field is a uniform flow
% u1 = unifrnd(-1,1); v1 =unifrnd(-1,1);
%
% %- Second flow: vortex cell flow
% Amp= (rand()<0.5)*unifrnd(-0.5,0.5); omiga= unifrnd(-0.5,0.5); phi= unifrnd(-pi,pi);
% u2 = @(x,y) Amp*sin(2*pi*omiga*x+phi).*cos(2*pi*omiga*y+phi);
% v2 = @(x,y)-Amp*cos(2*pi*omiga*x+phi).*sin(2*pi*omiga*y+phi);
%
% %- Third flow: a similar boundary flow (Mode-ratio bootstrapping method for PIV outlier correction)
% Amp3 = (rand()<0.5)*unifrnd(-0.5,0.5); x0 = unifrnd(-0.5,0.5); y0 = unifrnd(-0.5,0.5);
% u3 = @(x,y) Amp3*(x-x0).^2;
% v3 = @(x,y)-2*Amp3*(x-x0).*(y-y0);
%
% VecField(1,:) = u1+u2(Pos(1,:),Pos(2,:))+u3(Pos(1,:),Pos(2,:));
% VecField(2,:) = v1+v2(Pos(1,:),Pos(2,:))+v3(Pos(1,:),Pos(2,:));
%
% Vx = u1+u2(0.5,0.5) + u3(0.5,0.5);
% Vy = v1+v2(0.5,0.5) + v3(0.5,0.5);
%
% coef = maxV./max(abs([Vx,Vy])); % Take the maximal constraint into account
% VecField = coef*VecField; Vx = coef*Vx; Vy = coef*Vy;
end

% %% Add some noise to the particle images
% function I = AddNoise(Iin,xg,yg,opt)
% I = Iin + unifrnd(0,0.05)*randn(size(Iin));
% Ibackgrond = unifrnd(0,0.3)*exp(((xg-unifrnd(0,opt.Ni(1))).^2 +(yg-unifrnd(0,opt.Ni(2))).^2)/(-2*unifrnd(10,30)^2));
% I = max(I,Ibackgrond);
% end

function Iout = AddNoise(Iin,xg,yg,opt)
% Gaussian back ground for each pixel
Iout = Iin;
if rand<0.5 % Half probability
    gBackground = zeros(size(Iin),'single');
    gBackground = unifrnd(0,1/20)*rand(size(Iin))+ unifrnd(0,1/60);
    
    % Spatial illumination variance
    L = 1.4*max(opt.Ni);
    L0 = unifrnd(-0.5*L,1.5*L); sigma = unifrnd(L/10,L); angle = unifrnd(0,pi); Amp = unifrnd(0.0,0.25);
    sBackground = Amp.*exp(-(xg*sin(angle)+yg*cos(angle)-L0).^2./(2*sigma^2));
    
    Iout = max(max(Iin,gBackground),sBackground);
end
Iout = Iout + unifrnd(0,0.05)*randn(size(Iin));
end

%- Copy the field in opt to opts
function opts = argparse(opts, opt)
opts = opts;
fieldName = fieldnames(opt);
for i = 1:numel(fieldName)
    eval(cell2mat(['opts.',fieldName(i),' =', 'opt.',fieldName(i),';']))
end
end

%- Read the images from images and clip a patch
function  [I1,I2] = clipAndDefor(u,v,xg,yg,opts)
%- The path to store our original images
path = '../data/ImagesForDataset/';
Files = dir(fullfile(path,'*.bmp'));
LengthFiles = length(Files);
%- Read a random image
Img = single(imread(strcat(path,Files(randi(LengthFiles)).name)))/255;
[sx,sy] = size(Img);
%- Clip a random image patch
Idx = randi(sx-opts.Ni(1));Idy = randi(sy-opts.Ni(2));
Itemp = Img(Idx:(Idx+opts.Ni(2)-1),Idy:(Idy+opts.Ni(1)-1));

%- create the displacement field for the image
xn = (xg-1)/(opts.Ni(1)-1)-0.5;%- normalise the coordinate to [-0.5,0.5]
yn = (yg-1)/(opts.Ni(2)-1)-0.5;

Vx = u(xn,yn); Vy = v(xn,yn);
I1 = interp2(yg,xg,Itemp,yg+0.5*Vy,xg+0.5*Vx,'spline');
I2 = interp2(yg,xg,Itemp,yg-0.5*Vy,xg-0.5*Vx,'spline');
end
