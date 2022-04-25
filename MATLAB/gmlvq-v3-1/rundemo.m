% New for version 3.1 by Roland J. Veen, roland@rjv.at
%
function rundemo(number)
%RUNDEMO Run the specified demo (1-...)
%   number: number of the demo
%   1 - a seven class data set, single run with one prototype per class
%   See also run_single    
    fprintf('\nGMLVQ-OO 3.1 Demo\n--------------------------\n\n');
    demos = {
        'Demo I: a seven class data set, single run with one prototype per class'
        'Demo II: a simple two-class problem, 186-dim. feature vectors, 110 samples'
        'Demo III, IIIa: a difficult two-class problem, 32-dim. feature vectors, 98 samples'
        'Demo IV: reduced UCI segmentation data set (see demo III) validation run 10% rnd'
        'Demo V: a difficult two-class data set (see demo III) validation run 10% rnd'
        'Demo VI: a small subset of the simple two-class data (see demo II) leave one out'
        'Demo VII: UCI wine data set, classification according to high/low alcohol content'
    };
    if ~exist('number','var') % Missing parameter
        number = -1;
    end
    
    if strcmpi(number,'all'), number = 42; end 
    if strcmpi(number,'3a') || strcmpi(number,'iiia'), number = 31; end 
    
    if isa(number,'char') % Command syntax works also
        if contains(number, {'i','v','x'}, 'IgnoreCase', true)
            number = iRoman(number);
        else
            number = str2double(number);
        end
    end
    
    if (number > 0) && (number <= length(demos))
        fprintf('%s\n', demos{number});
    end
    switch number
        case 1
            load samplesData/uci-segmentation-sampled.mat fvec lbl;
            gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters(), 50);
            result = gmlvq.runSingle();
            result.plot();
        case 2            
            load samplesData/twoclass-simple.mat fvec lbl;
            gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters('rocClass', 1), 30);
            result = gmlvq.runSingle();
            plot(result); % Both calling styles work
        case 3
            load samplesData/twoclass-difficult.mat fvec lbl;
            gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters(), 50, [1 1 2]);
            result = gmlvq.runSingle();
            plot(result);
        case 31
            gmlvq = evalin('base','gmlvq');
            gmlvq.params.mu = 0.2;
            result = gmlvq.runSingle();
            plot(result);
        case 4
            load samplesData/uci-segmentation-sampled.mat fvec lbl;
            gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters(), 40);
            result = gmlvq.runValidation(10, 10);
            plot(result);
        case 5
            load samplesData/twoclass-difficult.mat fvec lbl;
            gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters(), 50, [1 1 2]);
            result = gmlvq.runValidation(10, 10);
            plot(result);
        case 6            
            load samplesData/twoclass-simple-small.mat fvec lbl;
            gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters(), 30);
            result = gmlvq.runL1O();
            plot(result);
        case 7
            fprintf('This is case 6 from 2.4, was missing from documentation v3.0, extrapolating...\n');
            % Wine-hi-lo added from 2.3...
            load samplesData/wine-hi-lo fvec lbl;
            gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters(), 50);
            result = gmlvq.runValidation(10, 10);
            plot(result);
        case 42
            n = 7;
            results = cell(1,n);
            for i = 1:n
                rundemo(i)
                results{i} = evalin('base', 'result');
            end
        otherwise
            if number == -1 
               fprintf('No number specified\n');
            else
                fprintf('Your input %i is not recognised.\n', number);
            end
            fprintf('Please choose one of the following:\n');
            for i = 1:length(demos)
                fprintf('%i: %s\n', i, demos{i});
            end
    end
    
    % Make data available in the workspace
    vars = { 'fvec' 'lbl' 'gmlvq' 'result' 'results' };
    for i = 1:length(vars)
        if exist(vars{i},'var'); assignin('base', vars{i}, eval(vars{i})); end
    end
end

function res = iRoman(number)
number = lower(number);
res = 0;
for i = 1:length(number)
    if number(i) == 'i' && i < length(number) && contains(number(i+1:length(number)), {'v','x'})
        res = res - 1;
    elseif number(i) == 'i'
        res = res + 1;
    elseif number(i) == 'v' && i < length(number) && contains(number(i+1:length(number)), {'x'})
        res = res - 5;
    elseif number(i) == 'v'
        res = res + 5;
    elseif number(i) == 'x'
        res = res + 10;
    end
end
end