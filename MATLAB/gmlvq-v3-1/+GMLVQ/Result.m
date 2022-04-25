classdef Result
    % Class that represents the result of a single algorithm run
    properties
        run GMLVQ.Run
        stepsizeMatrix (:,1) { mustBeNumeric }
        stepsizePrototypes (:,1) { mustBeNumeric }
    end
    
    methods
        % Plot implementation
        function plot(this, performances)
            if nargin < 2 % We have not specified a data set properly, plot the training data
                this.plot(this.run.trainingPerf);
                return;
            end
            
            close all;
            
            finalPerformance = performances(end);
            rocClass = num2str(this.run.gmlvq.params.rocClass);
            
            scrsz = get(0,'ScreenSize');
%            f1 = figure('Position',[1 scrsz(4)*8/10 scrsz(3)*4/10 scrsz(4)*8/10])   
            f1 = figure('Position',[1 scrsz(4)-(scrsz(4)*2/3) scrsz(3)/3 scrsz(4)*2/3]);   
            f1.Name = 'Overview';
            movegui(f1,'northwest'); %RJV FIX Y axis positioning
    
            % figure(1);                                       % learning curves 
            msize=15;   % size of symbols
            totl=this.run.nSteps+1;
            subplot(3,2,1);  % plot glvq cost fucntion vs. steps
            plot(1:totl,[performances.costFunction],':b.','MarkerSize',msize); hold on;
    
            title('glvq costs per example w/o penalty term ', 'FontName','LucidaSans', 'FontWeight','bold'); 
            xlabel('gradient steps');
            axis([0 totl -1 1]); axis 'auto y';
 
            subplot(3,2,2);  % plot total training error vs. steps
            plot(1:totl, [performances.totalError], ':r.', 'Markersize', msize); hold on;
   
            title('total training error', 'FontName','LucidaSans', 'FontWeight','bold'); 
            xlabel('gradient steps'); 
            axis([0 totl 0 1]); axis 'auto y';
   
            subplot(3,2,3); % plot the class-wise errors vs. steps
            %plot(1,performances(1).classWise,'.','MarkerSize',msize+10);
            %hold on; FIX RJV
            %plot(1,performances(1).classWise,'w.','MarkerSize',msize+10);
            p = plot(1:totl,[performances.classWise],':ko','MarkerSize',msize/3); % .
            % FIX RJV Make colors consistent
            pcol = GMLVQ.Util.consistentColorList();
            for i = 1:length(p)
%                 p(i).MarkerSize = 5;
                p(i).Color = pcol{i}; % Color
                p(i).MarkerEdgeColor = 'k';
                p(i).MarkerFaceColor = pcol{i}; % Color
            end
%             legend(num2str([1:this.run.nClasses]'),'Location','northoutside', ...
%                 'Orientation','horizontal'); %#ok<NBRAK> %,'NumColumns',2);
            legend(num2str([1:this.run.nClasses]'),'Location','best'); %#ok<NBRAK> %,'NumColumns',2);
  
            title('class-wise training errors', 'FontName','LucidaSans', 'FontWeight','bold'); 
            xlabel('gradient steps'); 
            axis([0 totl 0 1]);  axis 'auto y';
   
            subplot(3,2,4);   % plot AUC (ROC) vs. steps
            plot(1:totl, [performances.auroc], ':k.','MarkerSize',msize); 
            axis([ 0 totl min(0.9, min([performances.auroc])) 1.05 ]); 
            title(['AUC(ROC), class ',rocClass,' vs. all others'], 'FontName','LucidaSans', 'FontWeight','bold');
 
   
            subplot(3,2,5);  % plot glvq cost fucntion vs. steps
            plot([1:totl], this.stepsizePrototypes, ':b.','MarkerSize',msize); hold on;
            plot([1:totl], this.stepsizeMatrix, ':r.','Markersize',msize);
            title('stepsizes', 'FontName','LucidaSans', 'FontWeight','bold'); 
            xlabel('gradient steps');
            legend('prototype','relevances','Location','NorthEast');
            axis([0 totl -1 1]); axis 'auto y';
   
   
            f2 = figure(2);             % display the ROC curve of the final classifier 
            f2.Name = 'ROC-AUC';
            fprnpc = finalPerformance.fpr(finalPerformance.thresholds==0.5); % false positive of NPC
            tprnpc = finalPerformance.tpr(finalPerformance.thresholds==0.5); % true  positive of NPC
   
            plot(finalPerformance.fpr,finalPerformance.tpr,'-'); hold on;
            plot(fprnpc,tprnpc,'ko','MarkerSize',10,'MarkerFaceColor','b');
            axis square;
            xlabel('false positive rate');
            ylabel('true positive rate');
            title(['Training ROC, class ',rocClass,' (neg.) vs. all others (pos.)'], 'FontName','LucidaSans', 'FontWeight','bold'); 
            plot([0 1],[0 1],'k:'); 
            legend(['AUC = ',num2str(finalPerformance.auroc)],'NPC', 'Location','SouthEast'); %FIX RJV: Move after last plot
            hold off;
 

            f3 = figure(3);                   % visualize prototyeps and lambda matrix 
            f3.Name = 'Prototypes and Relevance Matrix';
            this.run.plot();
            movegui('northeast'); % RJV FIX Gui
            
            f4 = figure(4);                 % visualize data set in terms of projections on the
            f4.Name = 'Data Visualization'; % leading eigenvectors of Lambda
            this.run.visu_2d();
            title('2d-visualization of the data set', 'FontName','LucidaSans', 'FontWeight','bold');
            movegui('southeast'); % RJV FIX Gui
            f4.Position(2) = f4.Position(2) + 35; % taskbar
            
            % Addition: RJV show confusion matrix
            f5 = figure(5);
            f5.Name = 'Confusion Matrix';
            this.run.visu_conf();
            
            f5.Position(4) = round((scrsz(4)) / 2.5);
            f5.Position(3) = f5.Position(4);
            movegui(f5,'south')
            f5.Position(2) = f5.Position(2) + 35; % taskbar
        end
        
    end
end

