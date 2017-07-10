%% Experiments 3 (Particle concentration)
%- Investigate the pivdnn performance with different particle concentration in the synthetic uniform flow cases.
%- Adopt the root mean square (RMS) error
%- Yong Lee (leeyong@hust.edu.cn)
%- 2016-12-18


function exp3()
close all; clear;
%- Configuration
opts.MonteCarlo = 100; % 100
opts.show = false; % you can set it as true with debug model to see the intermediate results

opts.size = [128,128]; % Image size
opts.dt   = 2.2; %particle diameter
opts.inPlaneRatio = 0.999; % The ratio of the in-plane particles

%- Add methods pathes
addpath('../PIV_DNN/','../PIVCC/','../pivlab/','../OpticalFlow/');

Displacement = 0:0.1:10;
ParticleConcen  = [0.005, 0.01, 0.05, 0.1,0.5]; %particle concentration (particles per pixel)
BiasErr = zeros(3,numel(Displacement)); RMSErr = zeros(3,numel(Displacement));
for i = 1:numel(Displacement)
    for k = 1:numel(ParticleConcen)
        pivdnn_Res = [];
        for j = 1:opts.MonteCarlo
            [i,k,j]
            %- Generating the synthetic vector field and images
            ud = @(x,y)  Displacement(i) + 0*x;
            vd = @(x,y)  0 + 0*x;
            opts.u  = ud; opts.v = vd;
            opts.ppp = ParticleConcen(k);
            [x,y,u,v,I1,I2] = synSimu(opts);
            
            %- Process the images with pivdnn
            opt.display = false;opt.GPU = true; opt.Cudnn = true;
            [x_pivdnn,y_pivdnn,u_pivdnn,v_pivdnn] =pivdnn(I1,I2,opt);
            
            pivdnn_Res  = [pivdnn_Res; (u_pivdnn(:)-Displacement(i))];
        end
        BiasErr(k,i) = mean(pivdnn_Res(:));
        RMSErr(k,i)  = sqrt(mean(pivdnn_Res(:).^2));
    end
end
save exp3; % save the results

%- plot the results
for i = 1:numel(ParticleConcen),legendStr{i} = ['\fontsize{14}C_p = ',num2str(ParticleConcen(i)),' ppp'] ;end
Colours = linspecer(numel(ParticleConcen));
lineStyle = {':+'; ':'; '-.';'--'; '-';};

%- Mean bias error
H1 = figure; hold on; box on; xlim([min(Displacement), max(Displacement)]);%ylim([-0.01 1.01]);set(gca,'yscale','log');
for i = 1:numel(ParticleConcen)
    plot(Displacement,BiasErr(i,:),lineStyle{i},'LineWidth',2, ...
        'color',Colours(i,:),'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',Colours(i,:));
end
H11 = legend(legendStr);
xlabel('\fontsize{16}Particle image displacement \deltax [pixel]');
ylabel('\fontsize{16}Mean error(Bias) [pixel] ');set(gca,'fontsize',16); set(H1,'position',[300 100 600 350]);
%- Root Mean Square error
H2 = figure; hold on; box on; xlim([min(Displacement), max(Displacement)]);%ylim([-0.01 1.01]);set(gca,'yscale','log');
for i = 1:numel(ParticleConcen)
    plot(Displacement,RMSErr(i,:),lineStyle{i},'LineWidth',2, ...
        'color',Colours(i,:),'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',Colours(i,:));
end
H11 = legend(legendStr);
xlabel('\fontsize{16}Particle image displacement \deltax [pixel]');
ylabel('\fontsize{16}RMS error [pixel] ');set(gca,'fontsize',16); set(H2,'position',[400 200 600 350]);
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