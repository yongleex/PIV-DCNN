function [net, info] = cnn_pivdnn(varargin)
% Train the net with MatCovNet

%- Setup the MatConvNet package
run(fullfile(fileparts(mfilename('fullpath')),...
    '..','matconvnet', 'matlab', 'vl_setupnn.m')) ;

opts.batchNormalization = true ;
opts.network = [] ;
opts.networkType = 'simplenn' ;
[opts, varargin] = vl_argparse(opts, varargin) ;

sfx = opts.networkType ;
if opts.batchNormalization, sfx = [sfx '-bnorm'] ; end
opts.expDir = fullfile(vl_rootnn, 'data', ['dcnn-baseline-' sfx]) ;
[opts, varargin] = vl_argparse(opts, varargin) ;

opts.dataDir = fullfile(vl_rootnn, 'example', 'dcnn') ;
opts.imdbPath = fullfile('imdb.mat');
opts.train = struct() ;
opts.train.gpus =1; % 1 change to the first GPU and [] to CPU mode
opts.learningRate = 0.001;
opts.numEpochs    = 200;
opts.noiseLevel = 0.00;

opts = vl_argparse(opts, varargin) ;
global noise_std; noise_std = opts.noiseLevel;
if ~isfield(opts.train, 'gpus'), opts.train.gpus = []; end;

% --------------------------------------------------------------------
%                                                         Prepare data
% --------------------------------------------------------------------

if isempty(opts.network)
  net = cnn_pivdnn_init('batchNormalization', opts.batchNormalization, ...
    'networkType', opts.networkType,'learningRate',opts.learningRate,'numEpochs',opts.numEpochs);
else
  net = opts.network ;
  opts.network = [] ;
end

if exist(opts.imdbPath, 'file')
  imdb = load(opts.imdbPath) ;
else
  error('Please generate the training data at first!\n');
end

net.meta.classes.name = arrayfun(@(x)sprintf('%d',x),1:10,'UniformOutput',false) ;

% --------------------------------------------------------------------
%                                                                Train
% --------------------------------------------------------------------

switch opts.networkType
  case 'simplenn', trainfn = @cnn_train ;
  case 'dagnn', trainfn = @cnn_train_dag ;
end

[net, info] = trainfn(net, imdb, getBatch(opts), ...
  'expDir', opts.expDir, ...
  net.meta.trainOpts, ...
  'errorFunction','none',... 
  opts.train, ...
  'val', find(imdb.images.set == 3)) ;
end

% --------------------------------------------------------------------
function fn = getBatch(opts)
% --------------------------------------------------------------------
switch lower(opts.networkType)
  case 'simplenn'
    fn = @(x,y) getSimpleNNBatch(x,y) ;
  case 'dagnn'
    bopts = struct('numGpus', numel(opts.train.gpus)) ;
    fn = @(x,y) getDagNNBatch(bopts,x,y) ;
end
end

% --------------------------------------------------------------------
function [images, labels] = getSimpleNNBatch(imdb, batch)
% --------------------------------------------------------------------
images = imdb.images.data(:,:,:,batch) ;
labels = imdb.images.vector(1:2,batch) ;

% data augmentation by rotation and flip
k = randi(4);
images =  rot90(images,k);
labels =  [0,-1;1 0]^k*labels;

if rand > 0.5
    images = fliplr(images); 
    labels =  [1,0;0 -1]*labels;
end

%- Add noise
global noise_std;
if rand()< 0.5
    noiseStd = noise_std*rand();
    Inoise = imnoise(zeros(size(images),'single'),'gaussian',0,noiseStd);
    images = images + Inoise;
end

% --------------------------------------------------------------------
function inputs = getDagNNBatch(opts, imdb, batch)
% --------------------------------------------------------------------
images = imdb.images.data(:,:,:,batch) ;
labels = imdb.images.velocity(1:2,batch) ;

k = randi(4);
images =  rot90(images,k);
labels =  [0,-1;1 0]^k*labels;

if rand > 0.5
    images = fliplr(images); 
    labels =  [1,0;0 -1]*labels;
end

%- Add noise
%global noise_std;
if rand()< 0.5
    noiseStd = noise_std*rand();
    Inoise = imnoise(zeros(size(images),'single'),'gaussian',0,noiseStd);
    images = images + Inoise;
end


if opts.numGpus > 0
  images = gpuArray(images) ;
end
inputs = {'input', images, 'label', labels} ;
end
end

