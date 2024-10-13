% identify earthquakes

imageDir = 'brvideo/noise0001/';
inputPatchSize = [60 60 60 1];
outPatchSize = [20 20 20 1];
load('training_data_noise00004/checkpoints/checkpoint138120_epoch20.mat');

volReader = @(x) matRead(x);
voldsTest = imageDatastore(imageDir, ...
    'FileExtensions','.mat','ReadFcn',volReader);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin prediction

for id = 1:5
    disp(['Processing test volume ' num2str(id)]);
    
    vol{id} = read(voldsTest);
    volSize = size(vol{id},(1:3));
    volPadded = vol{id}(:,3:102,3:102);
    
    [heightPad,widthPad,depthPad,~] = size(volPadded);
    
    height = heightPad - (inputPatchSize(1)-outPatchSize(1));
    width = widthPad - (inputPatchSize(1)-outPatchSize(1));
    depth = depthPad - (inputPatchSize(1)-outPatchSize(1));
    
    tempSeg = zeros([height,width,depth]);
    
    % Overlap-tile strategy for segmentation of volumes.
    for k = 1:outPatchSize(3):depthPad-inputPatchSize(3)+1
        for j = 1:outPatchSize(2):widthPad-inputPatchSize(2)+1
            for i = 1:outPatchSize(1):heightPad-inputPatchSize(1)+1
                patch = volPadded( i:i+inputPatchSize(1)-1,...
                    j:j+inputPatchSize(2)-1,...
                    k:k+inputPatchSize(3)-1,:);
                patchSeg = predict(net,patch);
                tempSeg(i:i+outPatchSize(1)-1, ...
                    j:j+outPatchSize(2)-1, ...
                    k:k+outPatchSize(3)-1) = patchSeg;
            end
        end
	k
    end
    
    % Crop out the extra padded region.
    tempSeg = tempSeg(1:height,1:width,1:depth);

    % Save the predicted volume result.
    predictedLabels{id} = tempSeg;
end

save('brvideo/noise0001/predictions_noise0003/trained20epochs.mat', 'predictedLabels', '-v7.3');

