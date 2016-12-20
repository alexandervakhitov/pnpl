function [mnames, p] = xdraw_main(xs,yrange,method_list,field,ti,lx,ly, ws)
box('on');
hold('all');

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
    XTickLabs = {};
    for k = 1:length(xs)
        if (mod(k,2) == 0)
            XTickLabs {k} = num2str(xs(k));
        else
            XTickLabs {k} = ' ';
        end
    end
    set(gca, 'XTickLabel', XTickLabs);
    
%        );
end
ylim(yrange);
xlim(xs([1 end]));
set(gca,'xtick',xs);
set(gca,'FontSize',9);
title(ti,'FontSize',9, 'FontWeight', 'normal');
xlabel(lx,'FontSize',9);
ylabel(ly,'FontSize',9);


end