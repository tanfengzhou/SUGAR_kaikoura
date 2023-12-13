%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for AI prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startdate=1;
enddate=31;
starthour=0;
endhour=23;

imageDir = '12_3d/';
ai_threshold = 1.0;
%load('mat/eq_fewer00004_60/checkpoints/eqtrained3DUNetlinear10-Aug-2021-01-54-00-Epoch-20.mat'); checkpoint140440_epoch20.mat These two files are identical to EQwatcher_trained_model_fewer.mat   
load('EQwatcher_trained_model_fewer.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters from step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

win = 6;
startlat = -44.2;
startlon = 171.0;
latgrid = 0.0359983679;
longrid = 0.05004016838;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main program starts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputPatchSize = [60 60 60 1];
outPatchSize = [20 20 20 1];
valid = 20;

for date=startdate:enddate

if date<10
	day=strcat('0', num2str(date))
else
	day=num2str(date)
end

for hour=starthour:endhour

volLocTest = strcat(imageDir, day, num2str(hour),'.mat');
volReader = @(x) matRead(x);
voldsTest = imageDatastore(volLocTest, ...
    'FileExtensions','.mat','ReadFcn',volReader);

for id = 1:1
    disp(['Processing test volume ' day num2str(hour)]);
   
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


%save(strcat(imageDir, day, num2str(hour), '_ailoc.mat'),'predictedLabels');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AI prediction to map location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = predictedLabels{1};
rm = imregionalmax(E);
candidates = find(rm);

n=0;
for i = 1:length(candidates)
    if E(candidates(i)) > ai_threshold
        n=n+1;
        dadizhen(n)=candidates(i);
        score(n)=E(candidates(i));
    end
end

if n==0
	fclose(fopen(strcat(imageDir, day, num2str(hour), '_ailoc.csv'), 'w'));
	continue;
end

ss = size(E);
[steps, ys, xs] = ind2sub(ss, dadizhen);

n=1;

for i = 1:length(steps)
    x1 = steps(i)+20;
    x2 = ys(i)+22;
    x3 = xs(i)+22;
    br = vol{id}(x1,x2,x3);
    t = steps(i)-1;
    y = ys(i)-1;
    x = xs(i)-1;
    lat = startlat + (y+valid+2) * latgrid;
    lon = startlon + (x+valid+2) * longrid;
    time = (t+valid)/2 + win/4;
    catalog(n,:) = [time, lat, lon, score(n), br];
    n = n + 1;
end

catalog = sortrows(catalog);
writematrix(catalog, strcat(imageDir, day, num2str(hour), '_ailoc.csv'))

clear('dadizhen', 'score', 'catalog');

end

end


