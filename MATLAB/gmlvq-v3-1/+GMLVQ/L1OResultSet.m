classdef L1OResultSet < GMLVQ.ResultSet
    % Class that extends the default resultset with custom plotting behaviour
    methods
        function plot(this)
            close all;
            
            % Assign all variables
            mcftra = [this.averageRun.trainingPerf.costFunction];
            scftra = [this.averageRun.trainingPerfStd.costFunction];
            mcfval = [this.averageRun.validationPerf.costFunction];
            scfval = [this.averageRun.validationPerfStd.costFunction];
            
            mtetra = [this.averageRun.trainingPerf.totalError];
            stetra = [this.averageRun.trainingPerfStd.totalError];
            mteval = [this.averageRun.validationPerf.totalError];
            steval = [this.averageRun.validationPerfStd.totalError];
            
            mauctra = [this.averageRun.trainingPerf.auroc];
            sauctra = [this.averageRun.trainingPerfStd.auroc];
            maucval = [this.averageRun.validationPerf.auroc];
            saucval = [this.averageRun.validationPerfStd.auroc];
            
            mcwtra = [this.averageRun.trainingPerf.classWise]';
            scwtra = [this.averageRun.trainingPerfStd.classWise]';
            mcwval = [this.averageRun.validationPerf.classWise]';
            scwval = [this.averageRun.validationPerfStd.classWise]';
            
            totalsteps = this.averageRun.nSteps + 1;
            onlyat = floor(linspace(1, totalsteps, 10));
            
            rocClass = num2str(this.averageRun.gmlvq.params.rocClass);
            
            f1 = figure(1);
            msize = 15; % Size of symbols
            if totalsteps < 50; msize = 20; end
            
            subplot(2,2,1);   % total training error rate 
            % axis([0 totalsteps 0 1.2*max(mteval)]); 
            plot(1:totalsteps,mtetra,':.','MarkerSize',msize); 
            axis([1 totalsteps 0 1.05]); % axis 'auto y';
%             axis tight; axis 'auto y'; 
            hold on; box on;
            title('total training error rate',...
                'FontName','LucidaSans', 'FontWeight','bold'); 
            xlabel('gradient steps');
            errorbar(onlyat,mtetra(onlyat),stetra(onlyat)/sqrt(this.nRuns),...
                'co','MarkerSize',1); 
            legend('training','Location','Best'); %FIX Move legend
            hold off;

            subplot(2,2,2);  % total training AUC(ROC)
            plot(1:totalsteps,mauctra,':.','MarkerSize',msize); 
            axis([1 totalsteps 0.8 1.05]); % axis 'auto y';
            hold on; box on;
            % plot(1:totalsteps,maucval,'g.'); 
            title(['AUC(ROC) w.r.t. to class ',rocClass,' vs. all others'],...
                'FontName','LucidaSans', 'FontWeight','bold'); 
            errorbar(onlyat,mauctra(onlyat),sauctra(onlyat)/sqrt(this.nRuns),'b.'); 
            legend('training','Location','Best'); %FIX Move after errorbar
            hold off;

            subplot(2,2,3);   % class-wise training errors
            p = plot(1:totalsteps,mcwtra,':o','MarkerSize',msize/3); 
            % FIX RJV Make colors consistent
            pcol = GMLVQ.Util.consistentColorList();
            for i = 1:length(p)
               p(i).Color = pcol{i}; 
               p(i).MarkerEdgeColor = 'k';
               p(i).MarkerFaceColor = pcol{i}; % Color
            end            
            
            title('class-wise training errors',...
                'FontName','LucidaSans', 'FontWeight','bold');
            xlabel('gradient steps');
            legend(num2str([1:this.averageRun.nClasses]'),'Location','Best');
            axis tight; axis 'auto y'; 
            hold on; box on;
            hold off;

            subplot(2,2,4);   % cost function (training)
            plot(1:totalsteps,mcftra,':.','MarkerSize',msize); 
            title('cost fct. w/o penalty term (training)',...
                'FontName','LucidaSans', 'FontWeight','bold');
            xlabel('gradient steps');
            axis tight; axis 'auto y'; 
            hold on; box on;
            hold off;

            %  single l1O roc of final gmlvq systems after 
            %  totalsteps gradient steps

            f2 = figure(2); 
            fprs = this.averageRun.validationPerf(end).fpr;
            tprs = this.averageRun.validationPerf(end).tpr;
            thresh = this.averageRun.validationPerf(end).thresholds;
            plot(fprs,tprs,'b-','LineWidth',2);
            hold on;
            plot((fprs(thresh==0.5)),...
                (tprs(thresh==0.5)),'ko',...
            'MarkerSize',10,'MarkerFaceColor','g');
            plot([0 1],[0 1],'k:'); 
            legend(['AUC= ',num2str(-trapz((fprs),(tprs)))],... %RJV: Legend after last plot
                'NPC performance','Location','SouthEast');
            xlabel('false positive rate');
            ylabel('true positive rate'); 
            axis square; 
            title(['Leave-One-Out ROC (class ',rocClass,' vs. all others)'],...
                'FontName','LucidaSans', 'FontWeight','bold'); 
            hold off;
            
            f3 = figure(3);
            this.averageRun.plot();
            
            %ADD RJV Confusion Plot
            f4 = this.averageRun.visu_conf();
            
            f5 = figure;
            f5.Name = 'Data Visualization';
            this.averageRun.visu_2d();

            f6 = figure;
            f6.Name = 'Mean Confusion Matrix';
            subplot(1,2,1)
            this.averageRun.plotconf('training');
            subplot(1,2,2);
            this.averageRun.plotconf('validation');
            f6.Position(4) = round(f6.Position(3) / 2); % make height = width/2;
            
            f1.Name = 'Overview';
            f2.Name = 'ROC-AUC';
            f3.Name = 'Prototypes and Relevance Matrix';
            
            % Make presentable to the researcher (unstack the pancakes)
            f1.Position(4) = round(f1.Position(4) * 1.5); % slightly larger
            movegui(f1,'northwest');
            movegui(f2,'north');
            movegui(f3,'northeast');            
            movegui(f4,'south');
            f4.Position(2) = f4.Position(2) + 35; % taskbar
            movegui(f5,'southeast');
            f5.Position(2) = f5.Position(2) + 35; % taskbar

            movegui(f6,'southwest');
            f6.Position(2) = f6.Position(2) + 35; % taskbar            
        end
    end
end

