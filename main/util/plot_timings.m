function plot_timings
    load('lines_pts_number_timing');
    Xarg = 2*ns;
    for i = 1:2
        method_list(i).avg_t_acc1 = method_list(i).avg_t_1 + method_list(i).avg_t_2;
        method_list(i).avg_t_acc2 = method_list(i).avg_t_2 + method_list(i).avg_t_4;
    end
    

    w= 300;
    h= 300;
    XLabel = 'n_p+n_l';
    close all;    
    figure('color','w','position',[0,0,w,h]);
    yrange = [0 2e-2];    
    xdrawgraph(Xarg,yrange,method_list,'avg_t_acc1','Mean processing time',...
    XLabel,'Time, s.', ones(2,1));

    figure('color','w','position',[0,0,w,h]);
    yrange = [0 4e-2];    
    xdrawgraph(Xarg,yrange,method_list,'avg_t_acc2','Mean solving time',...
    XLabel,'Time, s.', ones(2,1));
end