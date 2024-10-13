clear;
load('brvideo/noise0001/predictions_noise0003/trained20epochs.mat');
complete = readmatrix('~/training_test_SUGAR/11clean/input_catalog.csv');

valid = 20;
%latnum = 60;
%lonnum = 60;
%latgrid = 0.036052350956181256;
%longrid = 0.044087335519604066;
latgrid = 0.0359983679;
longrid = 0.05004016838;

for hhh = 6:6

    hhh
    nnn=0;

    for thre = 0.1:0.1:3.0
        thre
        nnn=nnn+1;
        E = predictedLabels{5};
        rm = imregionalmax(E);
        candidates = find(rm);
    
        n=0;
        for i = 1:length(candidates)
            if E(candidates(i)) > thre
                n=n+1;
                dadizhen(n)=candidates(i);
                score(n)=E(candidates(i));
            end
        end
    
        ss = size(E);
        [steps, ys, xs] = ind2sub(ss, dadizhen);
    
        n=1;
    
        for i = 1:length(steps)
            t = steps(i)-1;
            y = ys(i)-1;
            x = xs(i)-1;
            lat = -44.2 + (y+valid+2) * latgrid;
            lon = 171.0 + (x+valid+2) * longrid;
            time = (t+valid)/2 + 1.5;
            catalog(n,:) = [time, lat, lon, score(n)];
            n = n + 1;
        end
    
        prediction = sortrows(catalog);
        %writematrix(prediction, strcat('brvideo/noise0001/predictions/30', num2str(hhh), '_ailoc_50epochs.csv'))

    
        timestart = hhh*3600+10;
        timeend = (hhh+1)*3600+60-10;         %-10;
    
        n = 1;
    
        for i=1:length(complete)
            if complete(i,1) > timestart && complete(i,1) < timeend
                truth(n,:) = complete(i,:);
                truth(n,1) = truth(n,1) - hhh*3600;
                n = n + 1;
            end
        end
    
        match = 0;
        found = ones(length(truth),1);
    
        n = 1;
    
        for i=1:length(prediction)
            for j=1:length(truth)
                if abs(prediction(i,1)-truth(j,1)) < 2 && distance(prediction(i,2), prediction(i,3), truth(j,2), truth(j,3)) < 0.2
                    match = match + 1;
                    found(j) = 0;
                    
                    pipei = j;
    
                    cha_t(n) = abs(prediction(i,1)-truth(j,1));
                    cha_epi(n) = distance(prediction(i,2), prediction(i,3), truth(j,2), truth(j,3));
    
                    if j<length(truth)
                        for k=j+1:length(truth)
                            if abs(prediction(i,1)-truth(k,1)) < 2 && distance(prediction(i,2), prediction(i,3), truth(k,2), truth(k,3)) < cha_epi(n)
                                cha_t(n) = abs(prediction(i,1)-truth(k,1));
                                cha_epi(n) = distance(prediction(i,2), prediction(i,3), truth(k,2), truth(k,3));
                                pipei = k;
                            end
                        end
                    end
                    
                    truth(pipei,1)=0;
                    n=n+1;
    
                    break;            
                end
            end
            if j==length(truth)
                prediction(i,:)
            end
        end
    
        missindex = find(found);
        missmag =  zeros(length(missindex),1);
        for i = 1:length(missindex)
            missmag(i) = truth(missindex(i), 5);
        end
    
        total_truth = length(truth)
        total_earthquakes = length(prediction)
    
        true_positive = match/length(prediction)
        recall = match/length(truth)
    
        error_t = mean(cha_t)
        error_epi = mean(cha_epi) * 111
        
        f1 = 2 / (1/true_positive + 1/recall)
        
        chart(nnn,1)=thre;
        chart(nnn,2)=total_earthquakes;
        chart(nnn,3)=true_positive;
        chart(nnn,4)=recall;
        chart(nnn,5)=f1;
    
        clear('dadizhen', 'score', 'catalog', 'truth', 'cha_t', 'cha_epi');
    end
    writematrix(chart, strcat('brvideo/noise0001/predictions_noise0003/11', num2str(hhh), '_trained_20epochs.csv'));

end

