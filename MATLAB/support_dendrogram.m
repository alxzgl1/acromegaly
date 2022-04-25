%-------------------------------------------------------------------------------
% Function: see, niak_visu_dendrogram.m
%-------------------------------------------------------------------------------
function support_dendrogram(X)

if X(1,1) <= X(end,1)
	flag_sim = false;
else
  flag_sim = true;
  X(:,1) = -X(:,1);
end

% size
[m, ~] = size(X);

% firstly sort out the order that the entities appear on the x-axis and stor in sc
sc = X(m, [2, 3]);                   
for i = 1:(m-1)
  j = m-i;
  l = length(sc);
  sci = find(sc == X(j, 4));        
  if isempty(sci)
    sc = [X(j, [2, 3]), sc(1:l)];
  else
    sc(sci) = [];
    sc = [sc(1:sci-1), X(j, [2, 3]), sc(sci:l-1)];
  end
end

rg = max(X(:,1)) - min(X(:,1));
ymin = min(X(:,1)) - (rg)/10;
ymax = max(X(:,1)) + (rg)/10;
if ymin == ymax
  ymax = 1;
  ymin = 0;
end

% set axis ranges
axis([0, m+2, ymin, ymax]);                 
hold on;

% handle flag
if flag_sim
  list_y = (ymin:(rg/10):ymax);
  set(gca, 'YTick', list_y);
  tYList = [];
  for q = 1:length(list_y), tYList{q} = sprintf('%1.3f', -list_y(q)); end
  set(gca, 'YTickLabel', tYList);
end

yc = ymin * ones(1, m+1);                     % Set initial yc as x-axis
ct = (1:m+1);                                 % Set initial x co-ord for plotting
i = 1;
while (i <= m) && (X(i,1) ~= Inf)
  y10 = yc(X(i, 2));
  y20 = yc(X(i, 3));
  y2 = X(i, 1);                               % y2 = level of join
  x1 = ct(find(ct(X(i, 2)) == sc));           % Find plotting co-ords for
  x2 = ct(find(ct(X(i, 3)) == sc));           % the entities joined                
  plot([x1, x1], [y10, y2], 'k');                  % Plot first vertical line
  plot([x2, x2], [y20, y2], 'k');                  % Plot second vertical line
  plot([x1, x2], [y2, y2], 'k');                   % Plot horizontal connecting line
  ct = [ct, (x1+x2+0.0001)/2];                %
  sc = [sc, (x1+x2+0.0001)/2];                % Update ct, sc and yc
  yc = [yc, y2];                              %
  i = i+1;
end

end %% end

%-------------------------------------------------------------------------------
