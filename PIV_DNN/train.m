clear; close all; rng('default');

ForceSynthesize = false; % Set true, means re- generating the training data, ignoring the existence of training data
DataSize = 100; % the number setting for each training dataset.

learningRate = [0.002,0.0014,0.0008,0.0005,0.0003]; % learning rate
numEpochs    = [  100,  125,   150,   175,   200];
% numEpochs    = [  25,   30,   40,    50,   60];

%- The path to restore the trainint dataset
path.F1   = '../data/NetF1/';
path.F2   = '../data/NetF2/';
path.F3   = '../data/NetF3/';
path.F4_1 = '../data/NetF4-1/';
path.F4_2 = '../data/NetF4-2/';
path.F4_3 = '../data/NetF4-3/';

%% Generate Network F1 training data
if (~ exist([path.F1,'imdb.mat'],'file')) || ForceSynthesize
    opt.mvr = [0,12];                       % The velocity component range.
    disp(['Generating the 1th (NetF1) random artificial PIV images'])
    [images, meta] = genPIVImgDB(DataSize,opt);
    save([path.F1,'imdb.mat'],'images','meta','-v7.3');
    clear images meta opt;
end

%% Generate Network F2 training data
if (~ exist([path.F2,'imdb.mat'],'file')) || ForceSynthesize
    opt.mvr = [0,5]; opt.Ni = [89,89];      % The velocity component range.
    disp(['Generating the 2th (NetF2) random artificial PIV images'])
    [images, meta] = genPIVImgDB(DataSize,opt);
    save([path.F2,'imdb.mat'],'images','meta','-v7.3');
    clear images meta opt;
end

%% Generate Network F3 training data
if (~ exist([path.F3,'imdb.mat'],'file')) || ForceSynthesize
    opt.mvr = [0,1.2]; opt.Ni = [77,77];     % The velocity component range.
    disp(['Generating the 3th (NetF3) random artificial PIV images'])
    [images, meta] = genPIVImgDB(DataSize,opt);
    save([path.F3,'imdb.mat'],'images','meta','-v7.3');
    clear images meta opt;
end

%% Generate Network F4_1 training data
if (~ exist([path.F4_1,'imdb.mat'],'file')) || ForceSynthesize
    opt.mvr = [0,0.1]; opt.Ni = [77,77];     % The velocity component range.
    disp(['Generating the 4th (NetF3_1) random artificial PIV images'])
    [images, meta] = genPIVImgDB_F4(DataSize,opt);
    save([path.F4_1,'imdb.mat'],'images','meta','-v7.3');
    clear images meta opt;
end

%% Generate Network F4_2 training data
if (~ exist([path.F4_2,'imdb.mat'],'file')) || ForceSynthesize
    opt.mvr = [0,0.1]; opt.Ni = [77,77];     % The velocity component range.
    disp(['Generating the 4th (NetF3_2) random artificial PIV images'])
    [images, meta] = genPIVImgDB_F4(DataSize,opt);
    save([path.F4_2,'imdb.mat'],'images','meta','-v7.3');
    clear images meta opt;
end

%% Generate Network F4_3 training data
if (~ exist([path.F4_3,'imdb.mat'],'file')) || ForceSynthesize
    opt.mvr = [0,0.1]; opt.Ni = [77,77];      % The velocity component range.
    disp(['Generating the 5th (NetF3_3) random artificial PIV images'])
    [images, meta] = genPIVImgDB_F4(DataSize,opt);
    save([path.F4_3,'imdb.mat'],'images','meta','-v7.3');
    clear images meta opt;
end

%% Train the networks with learning rate decay
for i = 1:numel(numEpochs);
    NetF1 = cnn_pivdnn('imdbPath',[path.F1,'imdb.mat'],'expDir',path.F1,'learningRate',learningRate(i),'numEpochs',numEpochs(i));
    NetF2 = cnn_pivdnn('imdbPath',[path.F2,'imdb.mat'],'expDir',path.F2,'learningRate',learningRate(i),'numEpochs',numEpochs(i));
    NetF3 = cnn_pivdnn('imdbPath',[path.F3,'imdb.mat'],'expDir',path.F3,'learningRate',learningRate(i),'numEpochs',numEpochs(i));
    NetF4_1 = cnn_pivdnn('imdbPath',[path.F4_1,'imdb.mat'],'expDir',path.F4_1,'learningRate',learningRate(i),'numEpochs',numEpochs(i));
    NetF4_2 = cnn_pivdnn('imdbPath',[path.F4_2,'imdb.mat'],'expDir',path.F4_2,'learningRate',learningRate(i),'numEpochs',numEpochs(i));
    NetF4_3 = cnn_pivdnn('imdbPath',[path.F4_3,'imdb.mat'],'expDir',path.F4_3,'learningRate',learningRate(i),'numEpochs',numEpochs(i));
    
    %% Remove the loss layer in the networks
    NetF1.layers(end)   = [];   NetF2.layers(end) = [];    NetF3.layers(end) = [];
    NetF4_1.layers(end) = []; NetF4_2.layers(end) = []; NetF4_3.layers(end) = [];
    %% save the networks
    save(['Nets',num2str(i),'.mat'],'NetF1','NetF2','NetF3','NetF4_1','NetF4_2','NetF4_3','-v7.3')
end
%% save the final networks
save('Nets.mat','NetF1','NetF2','NetF3','NetF4_1','NetF4_2','NetF4_3','-v7.3')

%% Test the training results
pivdnn();

