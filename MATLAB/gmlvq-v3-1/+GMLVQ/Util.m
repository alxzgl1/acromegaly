classdef Util
    %UTIL Utility functions for the GMLVQ
    %   Static class with
    properties (Constant)
        DeepLearningToolbox = 'Deep Learning Toolbox';
    end
    methods (Static)
        function colorList = consistentColorList()
           colorList = {'b';'r';'g';'c';'m';'y';'#FF9900';'#006600';'k'}; 
%            colorList = {'#0072BD';'#D95319';'#EDB120';'#7E2F8E';...	
%                '#77AC30';'#4DBEEE';'#A2142F';'y';'m';'g';'k'};        
        end
        
        function [targets, outputs] = deconfuse(confmat)
            % DECONFUSE reconstruct targets and outputs from a confusion
            % matrix, to be used for e.g. plotconfusion.
            % 
            % deconfuse(confmat), wehere confmat is an NxN confusion
            % matrix.
            % 
            % See also: PLOTCONFUSION.
            nclasses = length(confmat);
            targets = categorical(zeros(1,sum(confmat,'all')));
            outputs = targets;
            %disp(confmat);   
            % rebuild elements with their classification back here
            idx = 1;
            for i = 1:nclasses
                for j = 1:nclasses
                    n = confmat(i,j);
                    targets(idx:idx+n-1) = categorical(i-1);
                    outputs(idx:idx+n-1) = categorical(j-1);
                    idx = idx + n;
                end
            end
        end
        
        function props = propSet(props, field, value)
            if isnan(props.(field))
                props.(field) = value;
            else
                if props.(field) == value
                    fprintf('Notice: setting %s to the default value has no effect\n', field);
                else
                    fprintf('Notice: %s set to %g instead of recommended %g for current mode\n', field, props.(field), value);
                end
            end
        end
        
        function h = confplot(confmat, atitle)
            %   CONFPLOT requires Deep Learning Toolbox 
            %   make sure the active figure can receive the plot
            %   before calling the function
            if nargin < 2, atitle = ''; end
            
            if ~GMLVQ.Util.hasFunction('confusionchart')
                % Alternative display
                fprintf('Unable to find the function confusionchart, using fallback function.\n');
                fprintf('For best results use the Deep Learning Toolbox.\n');
                h = GMLVQ.Util.confplot_rjv(confmat, atitle);
            else
                h = GMLVQ.Util.confplot_cm(confmat, strtrim([atitle ' Confusion Matrix']));
%                 figure;
%                 GMLVQ.Util.confplot_rjv(confmat, atitle);
            end
        end
        
        function h = confplot_rjv(confmat, atitle, useparula)
            
            % Title is optional
            if nargin < 2
                atitle = '';
            end
            
            atitle = strtrim([atitle ' Confusion Matrix']);
            
            n = length(confmat);
            
            if nargin < 3
                useparula = false;
            end
            
            % Calculate percentages
            sums = sum(confmat,2);
            percs = confmat;
            percimg = confmat;
            for i = 1:n
                for j = 1:n
                    percs(i,j) = confmat(i,j) / sums(i);
                    if useparula || i == j
                        percimg(i,j) = percs(i,j);
                    else
                        percimg(i,j) = 1.01 + percs(i,j);
                    end
                end
            end
            
            if useparula
                cmap = parula(100);
                range = [0 1];
            else % Create custom color map
                % off diagonal
                funge = cool(200);
                % make zero even lighter so any error stands out
                funge(101,:) = [0.8025 0.7975 1]; % orig 0.5025    0.4975    1.0000
                % add summer colors for diagonal
                cmap = [flip(summer(101)); funge(101:200,:)];
                range = [0 2.01];
            end
            
            % Create base image
            colormap(cmap);
            im = imagesc(percimg, range);            
            im.AlphaData = 0.4;
            ax = gca;
            h = gcf;
%             h.Units = 'normalized';
            hold on
            
            for x = 1:n
                for y = 1:n
%                     plot([x x]+0.5,[0 n]+0.5, 'k');
%                     plot([0 n]+0.5,[y y]+0.5, 'k');
                    if(confmat(y,x) > 0)
                        if x == y
                            col = '#006400';
                        else
                            col = '#8B0000';
                        end
                        text(x,y-0.05, num2str(confmat(y,x)), ... %beware of transposition! y-0.15
                            'HorizontalAlignment','center', ...
                            'VerticalAlignment','middle', ...
                            'FontWeight','bold','FontSize',11);
                        text(x,y+0.25, GMLVQ.Util.numperc(percs(y,x)), ... %beware of transposition!
                            'HorizontalAlignment','center', ...
                            'VerticalAlignment','middle', ...
                            'FontWeight','normal','FontSize',9,'Color', col);
%                         text(x,y, num2str(confmat(y,x)), ... %beware of transposition! y-0.15
%                             'HorizontalAlignment','center', ...
%                             'VerticalAlignment','middle', ...
%                             'FontWeight','bold','FontSize',11);
                    end
                end
            end
            
            for i = 1:n
                plot([i i]+0.5,[0 n]+0.5, 'k');
                plot([0 n]+0.5,[i i]+0.5, 'k');                
            end
            
            axis square
            xlabel('Predicted class');
            ylabel('True class');
            xticks(1:n);
            xtickangle(45);
            yticks(1:n);
            ax.TickLength = [0 0];
            hold off
            title(atitle);
        end
        
        
        function str = numperc(f)
             str = num2str(f*100,'%.1f');
             if any(str == '.') % we have a decimal
                 k = length(str);
                 while str(k) == '0'
                     k = k -1;
                 end
                 if str(k) == '.', k = k - 1; end
                 str = [str(1:k) '%'];
             end
        end
        function h = confplot_cm(confmat, atitle)
            cm = confusionchart(confmat);
            cm.Title = atitle;
            cm.ColumnSummary = 'column-normalized';
            cm.RowSummary = 'row-normalized';
            
            cm.DiagonalColor = [0 0.7412 0.4471];
            % [188 230 196]/255
            
            h = cm.Parent;
            
        end
        
        function h = confplot_pc(confmat, title)
            %   CONFPLOT requires Deep Learning Toolbox 
            if nargin < 2, title = ''; end

            % For plotconfusion, we need to first deconfuse the matrix
            [targets, outputs] = GMLVQ.Util.deconfuse(confmat);
            % Now we switch targets and outputs to generate
            % the confusion matrix transposed, e.g. predicted class horz
            % and true class vertical
            h = plotconfusion(outputs, targets, title);
            % Plotconfusion starts numbering classes at zero, we want to
            % start at one. Solution: manually +1 the labels here.
            xt = xticklabels;
            yt = yticklabels;
            for i = 1:length(xt)
               xt{i} = char(xt{i}+1);
               yt{i} = char(yt{i}+1);               
            end
            xticklabels(xt);
            yticklabels(yt);
            ax = GMLVQ.Util.getFigureAxes(h);
            % Fix labels for transposed matrix
            ax.XLabel.String = 'Predicted class';
            ax.YLabel.String = 'True class';
            figure(h)            
        end
        
        function figureToSubplot(fSrc, fDest, axDest)
            ax = squeeze(findall(fSrc, 'type', 'axes'));
            dp = get(axDest,'position');
            axcp = copyobj(ax, fDest);
            set(axcp,'Position',dp);
            delete(axDest);
            close(fSrc);         
        end
        
        function cmcToSubplot(fSrc, fDest, axDest)
            ax = squeeze(findall(fSrc, 'type', 'ConfusionMatrixChart'));
            if isempty(ax) %fallback
                ax = squeeze(findall(fSrc, 'type', 'axes'));
            end
                
            dp = get(axDest,'position');
            axcp = copyobj(ax, fDest);
            set(axcp,'Position',dp);
            delete(axDest);
            close(fSrc);         
        end
                 
        function ax = getFigureAxes(theFigure)
            ax = squeeze(findall(theFigure, 'type', 'axes'));
        end
        
        function bool = hasToolbox(name)
            v = ver;
            bool = any(strcmp(cellstr(char(v.Name)), name));
        end
        
        function bool = hasFunction(name)
            bool = ~isempty(which(name));
        end
    end
end

