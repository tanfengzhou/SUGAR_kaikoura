imageDir = 'eq_fewer00004_60';

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

load('lgraph_linear_scaling');

%analyzeNetwork(lgraph)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% training options

options = trainingOptions('adam', ...
    'MaxEpochs',20, ...
    'InitialLearnRate',5e-4, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',5, ...
    'LearnRateDropFactor',0.95, ...
    'ValidationData',dsVal, ...
    'ValidationFrequency',2000, ...
    'Verbose',true, ...
    'CheckpointPath','eq_fewer00004_60/checkpoints', ...
    'MiniBatchSize',miniBatchSize);
    
%'Plots','training-progress', ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start training 

doTraining = true;
if doTraining
    modelDateTime = datestr(now,'dd-mmm-yyyy-HH-MM-SS');
    [net,info] = trainNetwork(dsTrain,lgraph,options);
    save([imageDir '/checkpoints/traininginfo-' modelDateTime '.mat'], 'info')
    save([imageDir '/checkpoints/eqtrained3DUNetlinear' modelDateTime '-Epoch-' num2str(options.MaxEpochs) '.mat'],'net');
    %save([imageDir '/checkpoints/eqtrained-' modelDateTime 'everything.mat'])
end




