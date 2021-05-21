%% Install and configure the system

function install( )
clear; close all;
%- Configuration based on your computer, this is a basic mode which can be
%work in simple CPU mode.
opt.enableGpu = true;
opt.enableCudnn = true;

if ispc() % windowns 7
    opt.cudaRoot  = 'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v8.0';
    opt.cudnnRoot = 'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/cuDNN_v6.0';
else % Linux system (Ubuntu 16.04)
    opt.cudaRoot  = '/usr/local/cuda-8.0';
    opt.cudnnRoot = 'local/cudnn5';
end

%- Download, Unzip and Install the MatConvNet if neccessary
if ~exist('matconvnet/matlab/vl_compilenn.m','file');
    %- Path information
    name = 'matconvnet-1.0-beta24'; ext = '.tar.gz';
    URL = ['http://www.vlfeat.org/matconvnet/download/',name,ext];
    fileName = [name,ext];
    
    disp('Downloading MatConvNet......')
    %  urlwrite(URL,fileName);    % download
    websave(fileName, URL);
    disp('Downloading MatConvNet Finished!')
    
    disp('Unziping MatConvNet......')
    untar(fileName, 'Temp');   % Unzip
    disp('Unziping MatConvNet Finished!')
    
    disp('Moving MatConvNet......')
    CP1_FLAG = copyfile(['Temp/',name],'matconvnet','f');  % Copy it to right dir
    disp('Moving MatConvNet Finished!')
    
    
    disp('Compiling MatConvNet......')
    addpath('matconvnet/matlab/');
    vl_compilenn('enableGpu', opt.enableGpu, ... % Compile the MatConvNet
        'cudaMethod', 'nvcc', ...
        'cudaRoot', opt.cudaRoot, ...
        'enableCudnn', opt.enableCudnn, ...
        'cudnnRoot', opt.cudnnRoot,...
        'Verbose', 1,...
        'EnableImreadJpeg',false);
    disp('Compiling MatConvNet Finished');
    
    rmdir('Temp','s'); % Remove the Temp folder and download file
    delete(fileName);
else
    disp('MatConvNet Files have been detected!');
end

% Copy some modified files to MatConvNet
CP2_FLAG = copyfile('PIV_DNN/ChangesForMatConvNet/vl_simplenn.m','matconvnet/matlab/simplenn','f');
CP3_FLAG = copyfile('PIV_DNN/ChangesForMatConvNet/euclideanloss.m','matconvnet/matlab','f');
CP4_FLAG = copyfile('PIV_DNN/ChangesForMatConvNet/smoothL1.m','matconvnet/matlab','f');

if CP2_FLAG && CP3_FLAG && CP4_FLAG
    disp('Modification of PIV-DNN has been installed successfully!');    
else
    disp('Modification of PIV-DNN has NOT been installed successfully!!!');    
end
disp('Runing an example of PIVnet!( < 30s)');
addpath('PIV_DNN');
pivdnn();
disp('Runing an example of PIVnet Finished!');
disp('PIVnet has been installed successfully!');
end

