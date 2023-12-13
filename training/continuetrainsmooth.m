imageDir = 'training_data_noise00004';

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
dsTrain = transform(dsVal,@(patchIn)augment_regression(patchIn,dataSource));

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

load('training_data_noise00004/checkpoints/net_checkpoint__138120__2023_06_26__03_55_02.mat','net')
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
    save(['training_data_noise00004/checkpoints/checkpoint138120_epoch20.mat'],'net');
end



