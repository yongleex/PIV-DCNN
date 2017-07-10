%% visualization of the conv filters
function exp7()
close all; clear;
addpath('../PIV_DNN/','../PIVCC/','../pivlab/','../OpticalFlow/');
load('Nets.mat')
%% layer 1;i = 1
weights = NetF1.layers{1}.weights{1};
normalizedWeights = (reshape(mapminmax(reshape(weights,[25,40])')',[5,5,40])+1)/2;
normalizedWeights = padarray(normalizedWeights,[1,1,0],0,'post');

viewNormalizedWeights = ...
    [normalizedWeights(:,:,1),normalizedWeights(:,:, 9),normalizedWeights(:,:,17),normalizedWeights(:,:,25),normalizedWeights(:,:,33);
     normalizedWeights(:,:,2),normalizedWeights(:,:,10),normalizedWeights(:,:,18),normalizedWeights(:,:,26),normalizedWeights(:,:,34);
     normalizedWeights(:,:,3),normalizedWeights(:,:,11),normalizedWeights(:,:,19),normalizedWeights(:,:,27),normalizedWeights(:,:,35);
     normalizedWeights(:,:,4),normalizedWeights(:,:,12),normalizedWeights(:,:,20),normalizedWeights(:,:,28),normalizedWeights(:,:,36);
     normalizedWeights(:,:,5),normalizedWeights(:,:,13),normalizedWeights(:,:,21),normalizedWeights(:,:,29),normalizedWeights(:,:,37);
     normalizedWeights(:,:,6),normalizedWeights(:,:,14),normalizedWeights(:,:,22),normalizedWeights(:,:,30),normalizedWeights(:,:,38);
     normalizedWeights(:,:,7),normalizedWeights(:,:,15),normalizedWeights(:,:,23),normalizedWeights(:,:,31),normalizedWeights(:,:,39);
     normalizedWeights(:,:,8),normalizedWeights(:,:,16),normalizedWeights(:,:,24),normalizedWeights(:,:,32),normalizedWeights(:,:,40);
    ];
% size(normalizedWeights)
% figure;hist(normalizedWeights(:))
% figure; 
% subplot(121); imshow(weights(:,:,2,20),[]);
% subplot(122); imshow(normalizedWeights(:,:,40),[]);

H = figure('name','The filters in 1st CONV layer');
colormap(jet);
pcolor(double(viewNormalizedWeights'));axis off;
colorbar('ytick',  [0,0.25,0.5,0.75,1],'FontSize',16);
set(H,'position',[200 400 800 420]);

%% layer 2;i = 4
weights = NetF1.layers{4}.weights{1};
size(weights);
normalizedWeights = (reshape(mapminmax(reshape(weights,[25,1000])')',[5,5,1000])+1)/2;
normalizedWeights = padarray(normalizedWeights,[1,1,0],0,'post');

viewNormalizedWeights = ...
    [normalizedWeights(:,:,1),normalizedWeights(:,:, 9),normalizedWeights(:,:,17),normalizedWeights(:,:,25),normalizedWeights(:,:,33);
     normalizedWeights(:,:,2),normalizedWeights(:,:,10),normalizedWeights(:,:,18),normalizedWeights(:,:,26),normalizedWeights(:,:,34);
     normalizedWeights(:,:,3),normalizedWeights(:,:,11),normalizedWeights(:,:,19),normalizedWeights(:,:,27),normalizedWeights(:,:,35);
     normalizedWeights(:,:,4),normalizedWeights(:,:,12),normalizedWeights(:,:,20),normalizedWeights(:,:,28),normalizedWeights(:,:,36);
     normalizedWeights(:,:,5),normalizedWeights(:,:,13),normalizedWeights(:,:,21),normalizedWeights(:,:,29),normalizedWeights(:,:,37);
     normalizedWeights(:,:,6),normalizedWeights(:,:,14),normalizedWeights(:,:,22),normalizedWeights(:,:,30),normalizedWeights(:,:,38);
     normalizedWeights(:,:,7),normalizedWeights(:,:,15),normalizedWeights(:,:,23),normalizedWeights(:,:,31),normalizedWeights(:,:,39);
     normalizedWeights(:,:,8),normalizedWeights(:,:,16),normalizedWeights(:,:,24),normalizedWeights(:,:,32),normalizedWeights(:,:,40);
    ];
% size(normalizedWeights)
% figure;hist(normalizedWeights(:))
% figure; 
% subplot(121); imshow(weights(:,:,2,20),[]);
% subplot(122); imshow(normalizedWeights(:,:,40),[]);

H = figure('name','The filters in 2nd CONV layer');
colormap(jet);
pcolor(double(viewNormalizedWeights'));axis off;
colorbar('ytick',  [0,0.25,0.5,0.75,1],'FontSize',16);
set(H,'position',[300 300 800 420]);


%% layer 3;i = 7
weights = NetF1.layers{7}.weights{1};
size(weights);
normalizedWeights = (reshape(mapminmax(reshape(weights,[16,5000])')',[4,4,5000])+1)/2;
normalizedWeights = padarray(normalizedWeights,[1,1,0],0,'post');

viewNormalizedWeights = ...
    [normalizedWeights(:,:,1),normalizedWeights(:,:, 9),normalizedWeights(:,:,17),normalizedWeights(:,:,25),normalizedWeights(:,:,33);
     normalizedWeights(:,:,2),normalizedWeights(:,:,10),normalizedWeights(:,:,18),normalizedWeights(:,:,26),normalizedWeights(:,:,34);
     normalizedWeights(:,:,3),normalizedWeights(:,:,11),normalizedWeights(:,:,19),normalizedWeights(:,:,27),normalizedWeights(:,:,35);
     normalizedWeights(:,:,4),normalizedWeights(:,:,12),normalizedWeights(:,:,20),normalizedWeights(:,:,28),normalizedWeights(:,:,36);
     normalizedWeights(:,:,5),normalizedWeights(:,:,13),normalizedWeights(:,:,21),normalizedWeights(:,:,29),normalizedWeights(:,:,37);
     normalizedWeights(:,:,6),normalizedWeights(:,:,14),normalizedWeights(:,:,22),normalizedWeights(:,:,30),normalizedWeights(:,:,38);
     normalizedWeights(:,:,7),normalizedWeights(:,:,15),normalizedWeights(:,:,23),normalizedWeights(:,:,31),normalizedWeights(:,:,39);
     normalizedWeights(:,:,8),normalizedWeights(:,:,16),normalizedWeights(:,:,24),normalizedWeights(:,:,32),normalizedWeights(:,:,40);
    ];
% size(normalizedWeights)
% figure;hist(normalizedWeights(:))
% figure; 
% subplot(121); imshow(weights(:,:,2,20),[]);
% subplot(122); imshow(normalizedWeights(:,:,40),[]);

H = figure('name','The filters in 3th CONV layer'); 
colormap(jet);
pcolor(double(viewNormalizedWeights'));axis off;
colorbar('ytick',  [0,0.25,0.5,0.75,1],'FontSize',16);
set(H,'position',[400 200 800 420]);
end