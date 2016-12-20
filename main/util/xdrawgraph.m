function xdrawgraph(xs,yrange,method_list,field,ti,lx,ly, ws, do_legend, w, h)
%the legend is at upper right in default
if (nargin < 9)
   do_legend = 0; 
end
box('on');
hold('all');

p= zeros(size(method_list));

for i= 1:length(method_list)
    mname = method_list(i).name;
    if (ws(i) ~= 1.0)
        mname = [mname ' (' num2str(ws(i)) ')'];
    end
    mnames{i} = mname;
    p(i)= plot(xs,ws(i)*method_list(i).(field),...
        'marker',method_list(i).marker,...
        'color',method_list(i).color,...
        'markerfacecolor',method_list(i).markerfacecolor,...
        'displayname',mname, ...
        'LineStyle', method_list(i).linestyle, ...
        'LineWidth',2,'MarkerSize',4);
    
%        );
end
ylim(yrange);
xlim(xs([1 end]));
set(gca,'xtick',xs);

% title(ti,'FontSize',12,'FontName','Arial');
% xlabel(lx,'FontSize',11);
% ylabel(ly,'FontSize',11);

title(ti,'FontSize',9, 'FontWeight', 'normal');
xlabel(lx,'FontSize',9);
ylabel(ly,'FontSize',9);


% XTickLabs = {};
% for i = 1:length(xs)
%     if (mod(i,2) == 0)
%         XTickLabs {i} = int2str(xs(i));
%     else
%         XTickLabs {i} = ' ';
%     end
% end
% set(gca, 'XTickLabel', XTickLabs);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(gcf, ['figures\' field '.pdf'], '-dpdf','-r0');

if (do_legend == 2)
   legend(p, 2); 
end
if (do_legend == 1)
    hold off;    
    fig2 = figure('position',[100,200,w,h]);
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);    
    legendflex(p, mnames, 'ref', fig2, 'box', 'off', 'padding', [5 10 50], 'FontSize', 18);
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, ['figures\legend.pdf'], '-dpdf', '-r0');
end


return
