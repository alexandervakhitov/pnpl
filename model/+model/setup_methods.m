function method_list = setup_methods(ch_names, trial_num, it_num, is_planar)
    % compared methods    
    A= zeros(it_num, 1);
    B= zeros(trial_num, it_num);
    met_type = {1,1,1,3,3,3,4,2,2,2,1};
    name= {'Mirzaei','RPnL', 'Pluecker', 'EPnP\_GN', 'EPnP\_GNN', 'OPnP', 'OPnP*', 'EPnPL', 'OPnPL', 'DLT', 'OPnP\_naive'};
    f= {   @mirzWrapper, @PNLWrapper, @plueckerWrapper, @EPnP_GNN, @EPnP_GNN, @OPnP, @OPnP, @EPnPLS_GN, @OPnPL, @DLT, @OPnP_E};
    marker= { '^', 'v', 's', 'o', 'x', 'x', 'none', '+','none', '<','v'};
    color= {'r','g','b','m','m','b','k','k','k','g','b'};
    markerfacecolor=  color;%{'r','g','n','m','n','n','r','r','g'};
    linestyle= {'-','--',':','-','-','-',':','-','-',':','-'};  
    min_pt_num = {0,0,0,6,6,3,3,6,3,6,3};
    min_ln_num = {6,6,10,0,0,0,3,6,3,6,3};
    
    name_planar = {'Mirzaei','RPnL', 'EPnP\_Planar\_GN', 'OPnP', 'OPnP*', 'EPnPL\_Planar', 'OPnPL', 'DLT\_Planar'};
    f_planar = {   @mirzWrapper, @PNLWrapper, @EPnP_planar_GN, @OPnP, @OPnP, @EPnPLS_Planar_GN, @OPnPL, @DLT_planar};
    met_type_planar = {1,1,3,3,4,2,2,2};
    marker_planar = { '^', 'v',  'o',  'x', 'none', '+','none', '<'};
    color_planar = {'r','g','m','b','k','k','k','g'};
    markerfacecolor_planar =  color_planar;%{'r','g','n','m','n','n','r','r','g'};
    linestyle_planar = {'-','--','-','-',':','-','-',':'};
    min_pt_num_planar = {0,0,0,0,3,6,3,6};
    min_ln_num_planar = {6,6,0,0,3,6,3,6};
    
    if (is_planar)
        f = f_planar;
        name = name_planar;
        met_type = met_type_planar;
        marker = marker_planar;
        color = color_planar;
        markerfacecolor = markerfacecolor_planar;
        linestyle = linestyle_planar;
        min_pt_num = min_pt_num_planar;
        min_ln_num = min_ln_num_planar;
    end
    
    

    flags = zeros(length(name), 1);
    inds = [];
    for i = 1:length(ch_names)
        ind = find(strcmp(name, ch_names{i}));
        inds = [inds; ind];
    end
    
    marker = marker(inds);
    color = color(inds);
    markerfacecolor = markerfacecolor(inds);
    linestyle = linestyle(inds);
    name = name(inds);
    f = f(inds);
    met_type = met_type(inds);
    min_pt_num = min_pt_num(inds);
    min_ln_num = min_ln_num(inds);
    
    method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A, 'mean_foc', A, 'mean_reproj', A, ...
        'med_r', A, 'med_t', A, 'med_foc', A, 'med_reproj', A, 'r', B, 't', B,...
        'foc', B, 'reproj', B, 'marker', marker, 'color', color, ...
        'markerfacecolor', markerfacecolor, 'linestyle', linestyle,'met_type',met_type, 'min_pt_num', min_pt_num, 'min_ln_num', min_ln_num);

end