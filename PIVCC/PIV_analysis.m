function [x,y,u,v] =  PIV_analysis(I1,I2,pass,opt1,opt2,opt3,opt4)
%- An implementation of iterative deformation of PIV estimation
%- Ref: Scarano, F. (2002). Iterative image deformation methods in PIV. Measurement Science and Technology, 13(1).
%- The idea of implementation is based on the realization of PIVlib
%- Yong Lee (leeyong@hust.edu.cn)
%- 2017-04-06

%- I1, I2: The input images, should be gray-value
%- pass: The pass number of the multi-pass method
%- optx: configuration of each number

%- Example :
%  opt1.x_start = 5;opt1.y_start = 5;opt1.x_step = 8;opt1.y_step = 8;opt1.x_win = 16; opt1.y_win = 16;
%  opt2 = opt1; opt3 = opt1; pass = 3
%  [x,y,u,v] =  PIV_analysis(I1,I2,pass,opt1,opt2,opt3)

if nargin < 1
    clear;close all
    I1 = imread('vp1a.tif');
    I2 = imread('vp1b.tif');
end
if ~exist('pass','var'), pass =3; end

I1 = double(I1); I2 = double(I2);


%- Default initialization of the algorithm
opts.x_start = 5; % seting it larger than 1 can avoid the boundary defect
opts.y_start = 5;
opts.x_step = 8;
opts.y_step = 8;
opts.x_win = 16;
opts.y_win = 16;

%- Adopt the usr settings
opt1.name = []; optss(1) = argparse(opts, opt1);% Get the usr defined configuration.
opt2.name = []; optss(2) = argparse(opts, opt2);
opt3.name = []; optss(3) = argparse(opts, opt3);
opt4.name = []; optss(4) = argparse(opts, opt4);


%- Some parameters and variable that will be used multiple times
epsi = 0.05*exp(-[0:0.9:8]);% The variation interval epsi is decreased for each iteration loop, according to an exponential decay
[Yq,Xq] = meshgrid(1:size(I1,2),1:size(I1,1));




%- Simple pre-processing of the particle images
I1 = I1 - mean(I1(:)); I1 = I1/std(I1(:));
I2 = I2 - mean(I2(:)); I2 = I2/std(I2(:));

%- Calculate the [u,v] using a traditional method for the first pass
u = 0*Xq; v = 0*Xq;
[u,v,~,x,y] = piv_analysis_with_deformation(I1,I2,u,v,Xq,Yq,optss(1));

% The main loop of iterative deformation
for i =2:pass
    %- deal with the outliers and using spline deformation
    warning('off') %#ok<*WNOFF>
    u = smoothn(u,'robust');
    v = smoothn(v,'robust');
    u_full=interp2(y,x,u,Yq,Xq,'spline');
    v_full=interp2(y,x,v,Yq,Xq,'spline');
    warning('on') %#ok<WNON>
    
    [u,v,~,x,y] = piv_analysis_with_deformation(I1,I2,u_full,v_full,Xq,Yq,optss(i));
    %     figure;quiver(v,u);% debug quiver
end

if nargin < 1
    figure;
    quiver(y,x,v,u); set(gca,'ydir','reverse')
    xlabel('Y & v direction');ylabel('X & u direction') ;
    title('The vector field');
end
end

%% Conduct image distoration and PIV analysis with normalized cross-correlation (NCC)
%- Only the second image is deformed with window offset and gradient compensation(3 steps)
%- 1. slicing/cutting the image with given deformation (Full vector field)
%- 2. calculating the integral shift
%- 3. calculating the subpixel shift
function [u_out,v_out,maxR,x,y] = piv_analysis_with_deformation(I1,I2,u_in,v_in,Xq,Yq,opts)
%- Get the slicing index of the first image
[x_siz,y_siz] = size(I1);
x_Index_start = opts.x_start:opts.x_step:(x_siz-opts.x_win);
y_Index_start = opts.y_start:opts.y_step:(y_siz-opts.y_win);
[y_Index_start,x_Index_start] = meshgrid(y_Index_start,x_Index_start); % All the grids left-upper index
[vector_size_x,vector_size_y] = size(x_Index_start);
Number_of_vetors = numel(x_Index_start);

[y_Index_win,x_Index_win] = meshgrid(0:(opts.y_win-1),0:(opts.x_win-1)); % The index for a window with 2 indices(x,y)

image_patches_index_x = repmat(x_Index_win, [1, 1, Number_of_vetors])+repmat(permute(x_Index_start(:),[2,3,1]), [opts.x_win, opts.y_win, 1]);
image_patches_index_y = repmat(y_Index_win, [1, 1, Number_of_vetors])+repmat(permute(y_Index_start(:),[2,3,1]), [opts.x_win, opts.y_win, 1]);

%- Image deformation
I1_deform = interp2(Yq,Xq,I1,Yq-0.5*v_in,Xq-0.5*u_in,'spline'); % spline method is better
I1_deform(isnan(I1_deform)) = 0;

I2_deform = interp2(Yq,Xq,I2,Yq+0.5*v_in,Xq+0.5*u_in,'spline');
I2_deform(isnan(I2_deform)) = 0;

%- Cut the image with interpolation, the spline interpolation method is
%highly recommended in PIVlab  instead of the 11*11 Gaussian weighted
%interpolation in Refs
image1_cuts = interp2(I1_deform,image_patches_index_y,image_patches_index_x,'spline');
image2_cuts = interp2(I2_deform,image_patches_index_y,image_patches_index_x,'spline');

% figure; imshow(image1_cuts(:,:,1),[])% debug show
% figure; imshow(image2_cuts(:,:,1),[])

%- Conduct the FFT-based cross correlation
image1_cuts_energy = sqrt(repmat(sum(sum(image1_cuts.^2,1),2),[opts.x_win, opts.y_win, 1]));
image2_cuts_energy = sqrt(repmat(sum(sum(image2_cuts.^2,1),2),[opts.x_win, opts.y_win, 1]));

R = fftshift(fftshift(real(ifft2(fft2(image1_cuts).*conj(fft2(image2_cuts)))),1),2)...
    ./(image1_cuts_energy.*image2_cuts_energy);% This is the correlation coefficience map
R(R<0.05) = 0.05; % set it larger than 0.05, for log(R) operation in subpixel estimation
% plot(R(:))% debug plot
mR = R + 0.001*rand(size(R)); % avoid multi-peaks
maxR = max(max(mR,[],1),[],2);

%- Find the integral shift with condition (R/maxR == 1)
R_max_one = mR./repmat(maxR,[opts.x_win, opts.y_win, 1]);
% plot(R_max_one(:))% debug plot
R_max_position = find(R_max_one == 1);

[u_index,v_index,vector_number] = ind2sub([opts.x_win, opts.y_win,Number_of_vetors],R_max_position);
vector_u_map = -(1:opts.x_win)' + round((opts.x_win+1)/2); % Mapping the maximal position in correlation map to the image shift
vector_v_map = -(1:opts.y_win)' + round((opts.y_win+1)/2);

u_out_integral = vector_u_map(u_index);
v_out_integral = vector_v_map(v_index);

%- Find the sub-pixel shift with  "Gaussian peak fit" method
%- Ref: Raffel, M., Willert, C. E., & Kompenhans, J. (2007). Particle image velocimetry: a practical guide. Springer Science & Business Media. Page 160
u_index(u_index<2) = 2;  u_index(u_index>opts.x_win-1) = opts.x_win-1; %
v_index(v_index<2) = 2;  v_index(v_index>opts.y_win-1) = opts.y_win-1;
% R_i_j = R(R_max_position); % R[i,j]
R_i_j_position = sub2ind([opts.x_win, opts.y_win,Number_of_vetors],u_index,v_index,vector_number);
R_i_j = R(R_i_j_position);% R[i,j]
R_im1_j_position = sub2ind([opts.x_win, opts.y_win,Number_of_vetors],u_index-1,v_index,vector_number);
R_im1_j = R(R_im1_j_position);% R[i-1,j]
R_ip1_j_position = sub2ind([opts.x_win, opts.y_win,Number_of_vetors],u_index+1,v_index,vector_number);
R_ip1_j = R(R_ip1_j_position);% R[i+1,j]
R_i_jm1_position = sub2ind([opts.x_win, opts.y_win,Number_of_vetors],u_index,v_index-1,vector_number);
R_i_jm1 = R(R_i_jm1_position);% R[i,j-1]
R_i_jp1_position = sub2ind([opts.x_win, opts.y_win,Number_of_vetors],u_index,v_index+1,vector_number);
R_i_jp1 = R(R_i_jp1_position);% R[i,j+1]

u_subpixel = -(log(R_im1_j)-log(R_ip1_j))./(2*log(R_im1_j)-4*log(R_i_j)+2*log(R_ip1_j));
v_subpixel = -(log(R_i_jm1)-log(R_i_jp1))./(2*log(R_i_jm1)-4*log(R_i_j)+2*log(R_i_jp1));

%- The relative velocity output is the combination of the integral and subpixel fractures
u_out = u_out_integral + u_subpixel;
v_out = v_out_integral + v_subpixel;

%- Provide the coordinate of the vectors and Matrix representation
x = x_Index_start + (opts.x_win-1)/2;
y = y_Index_start + (opts.y_win-1)/2;

u_out = reshape(u_out,[vector_size_x,vector_size_y]) + interp2(u_in,y,x); % add the deformation part of image distoration
v_out = reshape(v_out,[vector_size_x,vector_size_y]) + interp2(v_in,y,x);
maxR = reshape(maxR,[vector_size_x,vector_size_y]);

end

%- Copy the field (and value) in opt to opts
%- Yong Lee 2016.12
function opts_out = argparse(opts, opt)
opts_out = opts;
fieldName = fieldnames(opt);
for i = 1:numel(fieldName)
    eval(cell2mat(['opts_out.',fieldName(i),' =', 'opt.',fieldName(i),';']))
end
end
