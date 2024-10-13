fds = fileDatastore('training_data_noise00004/preprocessedDataset/labelsTr/*.mat', 'ReadFcn', @importdata);
fullFileNames = fds.Files;
numFiles = length(fullFileNames);

datas = fileDatastore('training_data_noise00004/preprocessedDataset/imagesTr/*.mat', 'ReadFcn', @importdata);
dataNames = datas.Files;

for k = 1 : numFiles
    if rem(k,1000) == 0
    	fprintf('Now reading file %s\n', fullFileNames{k});
    end
    load(fullFileNames{k});
    load(dataNames{k});

    if size(labelone) == [60 60 60] & size(dataone) == [60 60 60]
	    continue;
    end

    padlabel = zeros(60,60,60);
    if size(labelone) == [0 0 0] | size(dataone) == [0 0 0]
	    delete(fullFileNames{k});
	    fprintf('Now deleting file %s\n', fullFileNames{k});
	    delete(dataNames{k});
	    fprintf('Now deleting file %s\n', dataNames{k});
    else
    	padlabel(21:40,21:40,21:40) = labelone(:,:,:);
    	labelone = padlabel;
    
    	save(fullFileNames{k}, 'labelone');
    end

end

