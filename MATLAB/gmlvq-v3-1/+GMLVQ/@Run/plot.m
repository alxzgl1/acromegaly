% Visualizes a trained GMLVQ system in terms of prototypes and relevance matrix
function plot(this)

% define colors for up to seven classes
symbstrings = GMLVQ.Util.consistentColorList(); % ['b';'g';'r';'c';'m';'y';'k']; FIX Color consistency 

% geometry of figure showing prototypes
if this.nPrototypes<=6; rows=3;  else; rows=4;   end

col=ceil(this.nPrototypes/rows)+1;  % number of columns
posdiag=col;              % position of diagonal element bar plot

posiw=[];                 % positions of prototype bar plots
for i=1:rows
    posiw=[posiw,(1:col-1)+(i-1)*(col)];
end

for iw=1:this.nPrototypes                   %RJV TODO: Achieve color consitency between Here [@]
    % display prototypes as bar plots
    subplot(rows,col,posiw(iw)); 
    hold on;
    bar(this.prototypes(iw,:),'FaceColor',symbstrings{this.gmlvq.plbl(iw)}); 
    title(['prot. ',num2str(iw), ', class ',num2str(this.gmlvq.plbl(iw))], ... 
        'FontName','LucidaSans', 'FontWeight','bold'); 
    axis([0.3 this.trainingData.nDimensions+0.7 -max(max(abs(this.prototypes))) +max(max(abs(this.prototypes)))]); 
    grid on; 
    box on; 
end

% display diagonal matrix elements as bar plot
subplot(rows,col,posdiag); 
hold on;
bar(sort(eig(this.lambda),'descend')); title('eigenvalues of rel.mat.', ... 
    'FontName','LucidaSans', 'FontWeight','bold'); 
xlabel('eigenvalue index'); %FIX RJV 
grid on;  axis 'auto y'; box on;
axis([0.3 this.nDimensions+0.7 0 (0.01+max(diag(this.lambda)))]); 
hold off;

subplot(rows,col,posdiag+col);
hold on;
bar(diag(this.lambda)); title('rel. matrix, diag.', ... 
    'FontName','LucidaSans', 'FontWeight','bold'); 
xlabel('feature number'); 
grid on;  axis 'auto y'; box on;
axis([0.3 this.nDimensions+0.7 0 (0.01+max(diag(this.lambda)))]); 
hold off;

% display off-diagonal matrix elements as matrix
subplot(rows,col,3*col:col:col*rows);
lambdaoff = this.lambda.*(1-eye(this.nDimensions));     % zero diagonal 
imagesc(lambdaoff); colormap(summer(64)); % FIX Colormap number 
axis square; 
xlabel('off-diag. el.', ... 
    'FontName','LucidaSans', 'FontWeight','bold'); 
if col == 2 
    colorbar('Location','EastOutside');
else
    colorbar('Location','SouthOutside');
end
hold off;

end