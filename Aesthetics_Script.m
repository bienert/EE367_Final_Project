hAxes=findall(gcf,'type','axes');

if exist('hXlabel','var')
    set(hXlabel , ...
        'FontName'   , 'Times New Roman');

    set(hXlabel, ...
        'interpreter'   , 'tex', ...
        'FontSize'      , 13          );
end

if exist('hYlabel','var')
    set(hYlabel , ...
        'FontName'   , 'Times New Roman');

    set(hYlabel , ...
        'interpreter'   , 'tex', ...
        'FontSize'      , 13          );
end

if exist('hColorbar','var')
        set(hColorbar ,...
        'FontName'   , 'Times New Roman',...
        'FontSize'   , 10);
        hColorbarLabel=get(hColorbar,'YLabel');
        set(hColorbarLabel,'FontName','Times New Roman');
        set(hColorbarLabel,'FontSize',10);
end


if exist('hLegend','var')
    set(hLegend,'FontName','Times New Roman','FontSize'   , 10);
    set(hLegend , ...
    'interpreter'   , 'tex');
    hAxes=hAxes(hAxes~=hLegend);
    %set(hLegend,'Location','northeastoutside');
end

if exist('hTitle','var')
set( hTitle                    , ...
    'FontSize'   , 16          , ...
    'FontWeight' , 'normal'      );
set(hTitle , ...
    'FontName'   , 'Times New Roman');
set(hTitle , ...
    'interpreter'   , 'tex');
end

set( hAxes                       , ...
    'FontName'   , 'Times New Roman' );


set(hAxes, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1         );
% 
% Adjust linewidth
% hLines = get(hAxes,'children');
% for n=1:length(hLines);
%     if strcmpi(get(hLines(n),'type'),'line')
%         set(hLines(n),'LineWidth',0.75);
%     end
% end
