imageDir = 'training_data_noise000007_all3areas';

preprocessDataLoc = [imageDir filesep 'preprocessedDataset'];

volReader = @(x) matRead(x);
volLoc = fullfile(preprocessDataLoc,'imagesTr');
volds = imageDatastore(volLoc, ...
    'FileExtensions','.mat','ReadFcn',volReader);

lblLoc = fullfile(preprocessDataLoc,'labelsTr');
pxds = imageDatastore(lblLoc, ...
    'FileExtensions','.mat','ReadFcn',volReader);    
    
    
patchSize = [60 60 60];
patchPerImage = 1;
miniBatchSize = 25;
patchds = randomPatchExtractionDatastore(volds,pxds,patchSize, ...
    'PatchesPerImage',patchPerImage);
patchds.MiniBatchSize = miniBatchSize;

volLocVal = fullfile(preprocessDataLoc,'imagesVal');
voldsVal = imageDatastore(volLocVal, ...
    'FileExtensions','.mat','ReadFcn',volReader);

lblLocVal = fullfile(preprocessDataLoc,'labelsVal');
pxdsVal = imageDatastore(lblLocVal, ...
    'FileExtensions','.mat','ReadFcn',volReader);

dsVal = randomPatchExtractionDatastore(voldsVal,pxdsVal,patchSize, ...
    'PatchesPerImage',patchPerImage);
dsVal.MiniBatchSize = miniBatchSize;

dataSource = 'Training';
dsTrain = transform(patchds,@(patchIn)augment_regression(patchIn,dataSource));

dataSource = 'Validation';
dsVal = transform(dsVal,@(patchIn)augment_regression(patchIn,dataSource));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up 3-D Unet layers

inputPatchSize = [60 60 60 1];
encoderDepth = 2;
filters = 32;
outPatchSize = [20 20 20 1];

%analyzeNetwork(lgraph)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% training options

load('training_data_noise000007_all3areas/checkpoints/net_checkpoint__22195__2024_01_11__01_38_13.mat','net')
lgraph = layerGraph(net);

options = trainingOptions('adam', ...
    'MaxEpochs',1, ...
    'InitialLearnRate',0.00000000000000000000000000000000000000000000000001, ...
    'Verbose',true, ...
    'ValidationData',dsVal, ...
    'ValidationPatience',1, ...
    'MiniBatchSize',miniBatchSize);

%'Plots','training-progress', ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start training 

doTraining = true;
if doTraining
    modelDateTime = datestr(now,'dd-mmm-yyyy-HH-MM-SS');
    [net,info] = trainNetwork(dsTrain,lgraph,options);
    save(['training_data_noise000007_all3areas/checkpoints/checkpoint22195_epoch5.mat'],'net');
end

















