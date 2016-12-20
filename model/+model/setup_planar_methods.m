function method_list = setup_planar_methods(ch_names, trial_num, it_num)
    % compared methods    
    A= zeros(it_num, 1);
    B= zeros(trial_num, it_num);
    met_type = {1,1,1,3,3,3,4,2,2};
    name= {'Mirzaei','RPnL', 'EPnP\_Planar\_GN', 'OPnP', 'OPnP*', 'EPNPLS\_Planar\_GN', 'OPnPL'};
    f= {   @mirzWrapper, @PNLWrapper, @EPnP_planar_GN, @OPnP, @OPnP, @EPnPLS_Planar_GN, @OPnPL};
    marker= { '+', 'v', 's', 'd', 'o', 'p', '^', '>','<'};
    color= {'r','g','b','c','m','b','k','r','g'};
    markerfacecolor=  color;%{'r','g','n','m','n','n','r','r','g'};
    linestyle= {'-','--',':','-.','-','--',':','-.','-'};

    flags = zeros(length(name), 1);
    for i = 1:length(name)
        for j = 1:length(ch_names)
            if (strcmp(ch_names{j}, name{i}))
                flags(i) = 1;
            end
        end
    end
    inds = find(flags);
    marker = marker(inds);
    color = color(inds);
    markerfacecolor = markerfacecolor(inds);
    linestyle = linestyle(inds);
    name = name(inds);
    f = f(inds);
    met_type = met_type(inds);
    method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A, 'mean_foc', A, 'mean_reproj', A, ...
        'med_r', A, 'med_t', A, 'med_foc', A, 'med_reproj', A, 'r', B, 't', B,...
        'foc', B, 'reproj', B, 'marker', marker, 'color', color, ...
        'markerfacecolor', markerfacecolor, 'linestyle', linestyle,'met_type',met_type);

end