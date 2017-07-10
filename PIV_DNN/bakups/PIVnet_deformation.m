%% A PIV evaluation method using combinational CNNs (Convolutional Neural Networks)
function [x,y,u,v] =  PIVnet(I1,I2,opt)
%- u(down)?is the velocity along this x direction (down) in the image
%- v(right)?is the velocity along the horizontal y direction  in the image
if nargin == 0
    load ImgData;
end

% I1 = I1(1:100,:);
% I2 = I2(1:100,:);
% I1 = imread('cameraman.tif');
% I1 = I1(100:200,:);
%


%- Check the input data type
if isa(I1,'uint8'), I1 = single(I1)./255; I2 = single(I2)./255; end
if ~ isa(I1,'single'), error('Please provide the proper image format!\n'); end


%- Predine the default options
opts.x_bias = 5;
opts.y_bias = 5;
opts.x_stepSize = 16;
opts.y_stepSize = 16;


opts.multiPass = 3;
opts.GPU = false;
opts.Cudnn = false;
opts.batchSize = 100;


opts.stepSize = [16,16];
opts.minPosi  = [0,0]; % The interrogation windows w.r.t the image boundary
opts.padSize  = [20,20];


opts.epsilon=0.02;
opts.thresh=2;

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
net(3,1) = NetF3_1; net(3,2) = NetF3_2;
net(3,3) = NetF3_3; net(3,4) = NetF3_1;

opts.multiPass = size(net,1);

if  opts.GPU
    for i = 1:numel(net)
        net(i) = vl_simplenn_move(net(i),'gpu');
    end
end



%% - Initial operation
%- Get image size
[sxI,syI] = size(I1);         % image size

%- Get the interrogation window positions (left-up corner index (Ind_x,Ind_y))
x_Index_start = opts.x_bias: opts.x_stepSize: (sxI-64);
y_Index_start = opts.y_bias: opts.y_stepSize: (syI-64);
[y_Index_start,x_Index_start] = meshgrid(y_Index_start,x_Index_start); % All the grids left-upper index

%- Get the vector number
[vector_size_x,vector_size_y] = size(x_Index_start); % Vector field information
Number_of_vetors = numel(x_Index_start);

%- Initial velocities to zeros
u = zeros(vector_size_x,vector_size_y,'double'); v = u;

%- Pad the images to deal with the boundary cases
image1_roi = padarray(I1, opts.padSize, min(I1(:)),'both');
image2_roi = padarray(I2, opts.padSize, min(I1(:)),'both');
x_roi_Index_start = x_Index_start +  opts.padSize(1); % The roi index
y_roi_Index_start = y_Index_start +  opts.padSize(2);


%- Merge all the patches index
[y_Index_win,x_Index_win] = meshgrid(0:(64-1),0:(64-1));
image_roi_patches_index_x = repmat(x_Index_win, [1, 1, Number_of_vetors])+repmat(permute(x_roi_Index_start(:),[2,3,1]), [64, 64, 1]);
image_roi_patches_index_y = repmat(y_Index_win, [1, 1, Number_of_vetors])+repmat(permute(y_roi_Index_start(:),[2,3,1]), [64, 64, 1]);
linearInd_image_roi = sub2ind(size(image1_roi), image_roi_patches_index_x(:), image_roi_patches_index_y(:));
%- An example to get the image cuts with  linearInd_image_roi
% I1_crops_linear = image1_roi(linearInd_image_roi); I1_crops = reshape(I1_crops_linear,64,64,[]);

%- Predifine some constant matrix for speed-up the calculation
[yg_image,xg_image] = meshgrid(1:size(image1_roi,2),1:size(image1_roi,1));

%% - Main Net Loop
for i = 1:opts.multiPass
    %- Outlier detection and replacement
    [u,v] = NormMedFilter(u,v,opts.thresh,opts.epsilon);
    u = smoothn(u,'robust',1000);
    v = smoothn(v,'robust',1000);
    
    %- Get the dense vector field with interpolation, i.e., we can achieve the image deformation.
    y_temp = y_roi_Index_start+32;  y_temp = padarray(y_temp,[1,1],'replicate','both') ; y_temp(:,1) = 1;y_temp(:,end) =  size(xg_image,2) ;
    x_temp = x_roi_Index_start+32;  x_temp = padarray(x_temp,[1,1],'replicate','both') ; x_temp(1,:) = 1;x_temp(end,:) =  size(xg_image,1) ;
   
    u_dense = interp2(y_temp,x_temp,padarray(u,[1,1]),yg_image,xg_image,'spline');
    v_dense = interp2(y_temp,x_temp,padarray(v,[1,1]),yg_image,xg_image,'spline');
    figure; mesh(u_dense);
    
    u_dense(u_dense<min(u(:))) = min(u(:)); u_dense(u_dense>max(u(:))) = max(u(:)); v_dense(v_dense<min(v(:))) = min(v(:)); v_dense(v_dense>max(v(:))) = max(v(:));
%     u_dense = smoothn(u_dense,'robust',1000);
%     v_dense = smoothn(v_dense,'robust',1000);
    
    %- Image deformation operation
    image1_roi_deform = interp2(yg_image,xg_image,image1_roi,yg_image-0.5*v_dense,xg_image-0.5*u_dense,'spline');
    image2_roi_deform = interp2(yg_image,xg_image,image2_roi,yg_image+0.5*v_dense,xg_image+0.5*u_dense,'spline');
    
    %- Get the image patches according to current velocity
    image1_crops_linear = image1_roi_deform(linearInd_image_roi);
    image1_cut = reshape(image1_crops_linear,64,64,[]);
    image1_cut = image1_cut - repmat(mean(mean(image1_cut)),[64, 64, 1]); %- pre-processing the images, mean value substraction
    
    image2_crops_linear = image2_roi_deform(linearInd_image_roi);
    image2_cut = reshape(image2_crops_linear,64,64,[]);
    image2_cut = image2_cut - repmat(mean(mean(image2_cut)),[64, 64, 1]); %- pre-processing the images, mean value substraction
    
    %- Final image patches as a series of 2 channel images
    imgPatches =  cat(3,permute(image1_cut,[1 2 4 3]),permute(image2_cut,[1 2 4 3]));
    
    VectorTemp = zeros(1,1,2,Number_of_vetors);
    
    if  opts.GPU
        VectorTemp = gpuArray(VectorTemp);
    end
    
    %- Calculate the velocity residual using CNNs (gpu or cpu modes)
    %- DCNN_F1, DCNN_F2, DCNN_F3_1,
    for j=1:opts.batchSize:Number_of_vetors
        indxTemp = min(j+opts.batchSize-1,Number_of_vetors);
        if  opts.GPU
            Temp = vl_simplenn(net(i,1),gpuArray(imgPatches(:,:,:,j:indxTemp)),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
        else
            Temp = vl_simplenn(net(i,1),imgPatches(:,:,:,j:indxTemp),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
        end
        VectorTemp(:,:,:,j:indxTemp) = Temp(end).x;
        
    end
    Temp = permute(VectorTemp,[3 4 1 2]);
    uRes = reshape(Temp(1,:),[vector_size_x,vector_size_y]);
    vRes = reshape(Temp(2,:),[vector_size_x,vector_size_y]);
    
    if i == opts.multiPass
        %- DCNN_F3_1
        uRes1 = uRes; vRes1 = vRes;
        
        %- DCNN_F3_2
        for j=1:opts.batchSize:Number_of_vetors
            indxTemp = min(j+opts.batchSize-1,Number_of_vetors);
            if  opts.GPU
                Temp = vl_simplenn(net(i,2),gpuArray(imgPatches(:,:,:,j:indxTemp)),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            else
                Temp = vl_simplenn(net(i,2),imgPatches(:,:,:,j:indxTemp),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            end
            VectorTemp(:,:,:,j:indxTemp) = Temp(end).x;
        end
        Temp = permute(VectorTemp,[3 4 1 2]);
        uRes2 = reshape(Temp(1,:),[vector_size_x,vector_size_y]);
        vRes2 = reshape(Temp(2,:),[vector_size_x,vector_size_y]);
        
        %- DCNN_F3_3
        for j=1:opts.batchSize:Number_of_vetors
            indxTemp = min(j+opts.batchSize-1,Number_of_vetors);
            if  opts.GPU
                Temp = vl_simplenn(net(i,3),gpuArray(imgPatches(:,:,:,j:indxTemp)),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            else
                Temp = vl_simplenn(net(i,3),imgPatches(:,:,:,j:indxTemp),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            end
            VectorTemp(:,:,:,j:indxTemp) = Temp(end).x;
        end
        Temp = permute(VectorTemp,[3 4 1 2]);
        uRes3 = reshape(Temp(1,:),[vector_size_x,vector_size_y]);
        vRes3 = reshape(Temp(2,:),[vector_size_x,vector_size_y]);
        
        %- DCNN_F3_4
        for j=1:opts.batchSize:Number_of_vetors
            indxTemp = min(j+opts.batchSize-1,Number_of_vetors);
            if  opts.GPU
                Temp = vl_simplenn(net(i,4),gpuArray(imgPatches(:,:,:,j:indxTemp)),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            else
                Temp = vl_simplenn(net(i,4),imgPatches(:,:,:,j:indxTemp),[],[], 'mode','test','ConserveMemory',true,'CuDNN',opts.Cudnn);
            end
            VectorTemp(:,:,:,j:indxTemp) = Temp(end).x;
        end
        Temp = permute(VectorTemp,[3 4 1 2]);
        uRes4 = reshape(Temp(1,:),[vector_size_x,vector_size_y]);
        vRes4 = reshape(Temp(2,:),[vector_size_x,vector_size_y]);
        
        %- Get the average vector residual of these three results
        uRes = (uRes1+uRes2+uRes3+uRes4)/4;
        vRes = (vRes1+vRes2+vRes3+vRes4)/4;
        % uRes = (uRes1+uRes2+uRes3)/3;
        % vRes = (vRes1+vRes2+vRes3)/3;
    end
    
    %- Update the velocity by adding (the average of) the residuals
    u = u + uRes;
    v = v + vRes;
    if  opts.GPU, u = gather(u); v = gather(v);end
    
    x = x_Index_start +  32; % The roi index
    y = y_Index_start +  32;
    
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

h = fspecial('gaussian',5,4*sqrt(2));% alpha with equal weight
Temp1 = filter2(h,u.*OutlierLabels,'same');
Temp2 = filter2(h,OutlierLabels,'same')+eps;
uOut = Temp1./Temp2; uOut(OutlierLabels) = u(OutlierLabels);

Temp1 = filter2(h,v.*OutlierLabels,'same');
vOut = Temp1./Temp2; vOut(OutlierLabels) = v(OutlierLabels);
end

%- Copy the field in opt to opts
function opts = argparse(opts, opt)
opts = opts;
fieldName = fieldnames(opt);
for i = 1:numel(fieldName)
    eval(cell2mat(['opts.',fieldName(i),' =', 'opt.',fieldName(i),';']))
end
end
