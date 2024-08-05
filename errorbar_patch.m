function h = errorbar_patch(x,y,er,c,opto)

%ERRORBAR_PATCH    - errorbar by patch
%
% errorbar_patch(x,y,er,c)
%
%   input
%     - x:    
%     - y:    mean
%     - er:    
%     - color
%     - N:          start y value for plot
%     
%   output
%
if nargin < 4
    c = [0 0 1];
end

if size(x,1)~= 1
    x = x';
end

if size(y,1)~= 1
    y = y';
end

if size(x,1)~= 1
    er = er';
end

if nargin == 4
    opto = false;
end

if opto
    rectangle('Position', [0, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
    rectangle('Position', [1000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
    rectangle('Position', [2000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
    rectangle('Position', [3000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
end


X = [x fliplr(x)];
Y = [y+er fliplr(y-er)];
h1 = patch(X,Y,c,'edgecolor','none','FaceAlpha',0.2); hold on
% h1 = patch(X,Y,c,'edgecolor',c,'FaceAlpha',0.2); hold on
h2 = plot(x,y,'color',c,'LineWidth',1.5); hold on
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
if nargout>0, h = [h1 h2]; 
end

