function patchOut = augment_regression(patchIn,flag)

% Augment training data by randomly rotating and reflecting image patches.
% Do not augment validation data. For both training and validation data,
% crop the response to the network's output size. Return the image patches
% in a two-column table as required by the trainNetwork function for
% single-input networks.
%
% Copyright 2019 The MathWorks, Inc.
isValidationData = strcmp(flag,'validation');

inpVol = cell(size(patchIn,1),1);
inpResponse = cell(size(patchIn,1),1);

yes=0;
no=0;
n=1;

% 5 augmentations: nil,rot90,fliplr,flipud,rot90(fliplr)
fliprot = @(x) rot90(fliplr(x));
augType = {@rot90,@fliplr,@flipud,fliprot};

for id=1:size(patchIn,1) 
    rndIdx = randi(8,1);
    tmpImg =  patchIn.InputImage{id};
    tmpResp = patchIn.ResponseImage{id};
    if rndIdx > 4 || isValidationData
        out =  tmpImg;
        respOut = tmpResp;
    else
        out =  augType{rndIdx}(tmpImg);
        respOut = augType{rndIdx}(tmpResp);
    end
    % Crop the response to to the network's output.
    respFinal=respOut(21:end-20,21:end-20,21:end-20,:);
    if max(respFinal,[],'all') >= 0
        inpVol{n,1}= out;
        inpResponse{n,1}=respFinal;
        yes = yes + 1;
        n = n + 1;
    else
        if no <= 100
            inpVol{n,1}= out;
            inpResponse{n,1}=respFinal;
            no = no + 1;
            n = n + 1;
        end
    end
end
patchOut = table(inpVol,inpResponse);










