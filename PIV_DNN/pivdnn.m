%% A PIV evaluation method using deep neural networks (DNN) (Convolutional Neural Networks)
function [x,y,u,v] =  pivdnn(I1,I2,opt)
%- u(down)?is the velocity along this x direction (down) in the image
%- v(right)?is the velocity along the horizontal y direction  in the image
if nargin == 0
    load ImgData;
end
% I1 = I1(1:100,:);
% I2 = I2(1:100,:);

%- Check the input data type
if isa(I1,'uint8'), I1 = single(I1)./255; I2 = single(I2)./255; end
if ~ isa(I1,'single'), error('Please provide the proper image format!\n'); end

I1 = I1- mean(I1(:)); % mean value substraction
I2 = I2- mean(I2(:));

%- Predine the default options
opts.stepSize = [16,16];
opts.minPosi  = [0,0]; % The interrogation windows w.r.t the image boundary
opts.padSize  = [20,20];
opts.multiPass = 3;
opts.GPU = false;
opts.Cudnn = false;
opts.batchSize = 100;

opts.epsilon = 0.02;
opts.thresh  = 2;

opts.display = true;

opt.name = [];
opts = argparse(opts, opt);% Get the usr defined configuration.

%- Setup the MatConvNet package
run(fullfile(fileparts(mfilename('fullpath')),...
    '..','matconvnet', 'matlab', 'vl_setupnn.m')) ;
%- load the net
load('Nets.mat');

net(1,1) = NetF1;
net(2,1) = NetF2;
net(3,1) = NetF3;
net(4,1) = NetF4_1; net(4,2) = NetF4_2; net(4,3) = NetF4_3;
opts.multiPass = size(net,1);

if  opts.GPU
    for i = 1:numel(net)
        net(i) = vl_simplenn_move(net(i),'gpu');
    end
end

%% - Initial operation
%- Get the basic information about the images
[sxI,syI] = size(I1);         % image size
sizV = floor(([sxI-64-2*opts.minPosi(1),syI-64-2*opts.minPosi(2)])./opts.stepSize)+1;
sxV = sizV(1); syV = sizV(2); % size/number of vectors

u = zeros(sxV,syV,'single');  % Initial velocities
v = zeros(sxV,syV,'single');

%- Pad the images to deal with the boundary cases
image1_roi = padarray(I1, opts.padSize, min(I1(:)),'both');
image2_roi = padarray(I2, opts.padSize, min(I1(:)),'both');
[y_Index_Img_roi,x_Index_Img_roi] = meshgrid(1:size(image1_roi,2),1:size(image1_roi,1));

%- The pathches position limits
minix = opts.minPosi(1) + opts.padSize(1);
miniy = opts.minPosi(2) + opts.padSize(2);
maxix = minix+ opts.stepSize(1) * (sxV-1);
maxiy = miniy+ opts.stepSize(2) * (syV-1);

%- The start position of patches
start_x = minix:opts.stepSize(1):maxix;
start_y = miniy:opts.stepSize(2):maxiy;
[start_y,start_x] = meshgrid(start_y,start_x);
x = start_x-opts.padSize(1) + 32.5;
y = start_y-opts.padSize(2) + 32.5;
start_x = permute(start_x(:), [2 3 1]); % denote the left-up conner
start_y = permute(start_y(:), [2 3 1]); % denote the left-up conner

%- Get the initial index of each patches
x_index = repmat(start_x,[64 64 1]) + repmat(meshgrid(0:63)',[1,1,numel(start_x)]);
y_index = repmat(start_y,[64 64 1]) + repmat(meshgrid(0:63) ,[1,1,numel(start_x)]);

%% - Main Net Loop
for i = 1:opts.multiPass
    %- Outlier detection and replacement
    if i<=3
    [u,v] = NormMedFilter(u,v,opts.thresh,opts.epsilon);
    end
    
    u_p = repmat(permute(u(:), [2 3 1]), [64, 64, 1]);
    v_p = repmat(permute(v(:), [2 3 1]), [64, 64, 1]);
    
    %- Update the index of patches
    x_index1 = double(x_index - 0.5*u_p);
    y_index1 = double(y_index - 0.5*v_p);
    x_index2 = double(x_index + 0.5*u_p);
    y_index2 = double(y_index + 0.5*v_p);
    
    %- Get the patches after window offset
    image1_cut = interp2(y_Index_Img_roi,x_Index_Img_roi,image1_roi,y_index1,x_index1,'spline'); % get the patches in the first image
    image2_cut = interp2(y_Index_Img_roi,x_Index_Img_roi,image2_roi,y_index2,x_index2,'spline'); % get the patches in the second image
    imgPatches =  cat(3,permute(image1_cut,[1 2 4 3]),permute(image2_cut,[1 2 4 3]));
    
    VectorTemp = zeros(1,1,2,prod(sizV));
    
    if  opts.GPU
        VectorTemp = gpuArray(VectorTemp);
    end
    
    %- Calculate the velocity residual using CNNs (gpu or cpu modes)
    %- DCNN_F1, DCNN_F2, DCNN_F3_1,
    for j=1:opts.batchSize:prod(sizV)
        indxTemp = min(j+opts.batchSize-1,prod(sizV));
        if  opts.GPU
            Temp = vl_simplenn(net(i,1),gpuArray(imgPatches(:,:,:,j:indxTemp)),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
        else
            Temp = vl_simplenn(net(i,1),imgPatches(:,:,:,j:indxTemp),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
        end
        VectorTemp(:,:,:,j:indxTemp) = Temp(end).x;
        
    end
    Temp = permute(VectorTemp,[3 4 1 2]);
    uRes = reshape(Temp(1,:),sizV);
    vRes = reshape(Temp(2,:),sizV);
    
    if i == opts.multiPass
        %- DCNN_F3_1
        uRes1 = uRes; vRes1 = vRes;
        
        %- DCNN_F3_2
        for j=1:opts.batchSize:prod(sizV)
            indxTemp = min(j+opts.batchSize-1,prod(sizV));
            if  opts.GPU
                Temp = vl_simplenn(net(i,2),gpuArray(imgPatches(:,:,:,j:indxTemp)),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            else
                Temp = vl_simplenn(net(i,2),imgPatches(:,:,:,j:indxTemp),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            end
            VectorTemp(:,:,:,j:indxTemp) = Temp(end).x;
        end
        Temp = permute(VectorTemp,[3 4 1 2]);
        uRes2 = reshape(Temp(1,:),sizV);
        vRes2 = reshape(Temp(2,:),sizV);
        
        %- DCNN_F3_3
        for j=1:opts.batchSize:prod(sizV)
            indxTemp = min(j+opts.batchSize-1,prod(sizV));
            if  opts.GPU
                Temp = vl_simplenn(net(i,3),gpuArray(imgPatches(:,:,:,j:indxTemp)),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            else
                Temp = vl_simplenn(net(i,3),imgPatches(:,:,:,j:indxTemp),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            end
            VectorTemp(:,:,:,j:indxTemp) = Temp(end).x;
        end
        Temp = permute(VectorTemp,[3 4 1 2]);
        uRes3 = reshape(Temp(1,:),sizV);
        vRes3 = reshape(Temp(2,:),sizV);
        
        uRes = (uRes1+uRes2+uRes3)/3;
        vRes = (vRes1+vRes2+vRes3)/3;
    end
    
    %- Update the velocity by adding (the average of) the residuals
    u = u + uRes;
    v = v + vRes;
    if  opts.GPU, u = gather(u); v = gather(v);end
    
    %- Display the vector field for debug, make the coordinate consistency
    if opts.display, figure();quiver(y,-x,v,-u); title(['the ',int2str(i),'-th level vector prediction']);end
end


end


%- Normalized median filter for outlier removal
%- The detected outliers are replaced with the weighted Gaussian smooth value
%- Ref: Westerweel, J. & Scarano, F. Exp Fluids (2005) 39: 1096. doi:10.1007/s00348-005-0016-6
function [uOut,vOut] =  NormMedFilter(u,v,Thr,eps)

[J,I]=size(u);

NormFluct=zeros(J,I,2);

%- Deal with the boundary elements
b=2;
J =  J+2*b; I = I + 2*b;
U = padarray(u,[b b],'symmetric','both');
V = padarray(v,[b b],'symmetric','both');


for c=1:2
    if c==1;VelComp=U;else VelComp=V;end;
    
    for i=1+b:I-b
        for j=1+b:J-b
            Neigh = VelComp(j-b:j+b,i-b:i+b);
            NeighCol=Neigh(:);
            NeighCol2=[NeighCol(1:(2*b+1)*b+b);NeighCol((2*b+1)*b+b+2:end)];
            MedianV=median(NeighCol2);
            Fluct=VelComp(j,i)-MedianV;
            Res=NeighCol2-MedianV;
            MedianRes=median(abs(Res));
            NormFluct(j-b,i-b,c)=abs(Fluct/(MedianRes+eps));
        end;
    end;
end;

OutlierLabels =sqrt(NormFluct(:,:,1).^2+NormFluct(:,:,2).^2)<Thr;

h = fspecial('gaussian',7,4*sqrt(2));% alpha with equal weight
Temp1 = filter2(h,u.*OutlierLabels,'same');
Temp2 = filter2(h,OutlierLabels,'same')+eps;
uOut = Temp1./Temp2; %uOut(OutlierLabels) = u(OutlierLabels);

Temp1 = filter2(h,v.*OutlierLabels,'same');
vOut = Temp1./Temp2; %vOut(OutlierLabels) = v(OutlierLabels);


end

%- Copy the field in opt to opts
function opts = argparse(opts, opt)
opts = opts;
fieldName = fieldnames(opt);
for i = 1:numel(fieldName)
    eval(cell2mat(['opts.',fieldName(i),' =', 'opt.',fieldName(i),';']))
end
end
