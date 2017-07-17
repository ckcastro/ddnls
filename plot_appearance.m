function plot_appearance(gcf,filename)
% usage plot_appearance(gcf)
%      Input -     gcf : figure handle
%              filename: string of the filename
     
% Tighthen white space from figure.
% Save pdf files. 

fontsize_value = 25;

set(0,'defaultlinelinewidth',2) % line thickness everywhere
set(0,'defaultaxesfontsize',fontsize_value)

set(0,'DefaultTextFontSize',fontsize_value); % Global text fontsize
set(gca, 'FontSize', fontsize_value); % Axis thicks fontsize


% global legend fontsize
%set(0,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual')


set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);

filename = strcat(pwd,'/',filename);
print('-dpdf',filename);
%print('-depsc',filename)
%print('-djpeg',filename)

%savefig(filename)


end

