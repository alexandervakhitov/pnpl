function plot_lines_number
load('planar_lines_number.mat');
% 
% num = 100;
% method_names = {'Mirzaei','RPnL', 'Pluecker', 'EPnP\_GNN', 'OPnP', 'OPnP*', 'EPNPLS\_GN', 'OPnPL'};
% method_list = model.setup_methods(method_names, num, length(ns));

close all;
yrange= [0 5];

i= 0; w= 300; h= 300;

XLabel = 'Lines and points number';
Xarg = 2*ns;
ws = ones(length(method_list), 1);
yrange = [0 1.5];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(Xarg,yrange,method_list,'med_r','Median Rotation',...
    XLabel,'Rotation Error (degrees)', ws, 1, w, h);

yrange = [0 7];
ws = ones(length(method_list), 1);
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(Xarg,yrange,method_list,'mean_r','Mean Rotation',...
    XLabel,'Rotation Error (degrees)',ws, 2);

yrange= [0 2];
ws = ones(length(method_list), 1);
ws(3) = 0.25;
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(Xarg,yrange,method_list,'med_t','Median Translation',...
    XLabel,'Translation Error (%)',ws, 2);

yrange= [0 5];
ws = ones(length(method_list), 1);
ws(1) = 0.5;
ws(3) = 0.5;
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(Xarg,yrange,method_list,'mean_t','Mean Translation',...
    XLabel,'Translation Error (%)',ws, 2);
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(Xarg,yrange,method_list,'med_reproj','Reprojection',...
%     XLabel,'Reprojection Error (pixels)',2);
% 
% figure('color','w','position',[w*(i-4),200,w,h]);i=i+1;
% xdrawgraph(Xarg,yrange,method_list,'med_reproj_lines','Reprojection For Lines',...
%     XLabel,'Reprojection Error (pixels)',2);
% 
% figure('color','w','position',[w*(i-5),200,w,h]);i=i+1;
% xdrawgraph(Xarg,yrange,method_list,'med_reproj_pts_lines','Reprojection For Both',...
%     XLabel,'Reprojection Error (pixels)',2);
% 
yrange = [0 20];
ws = ones(length(method_list), 1);
ws(1) = 0.12;
ws(3) = 0.3;
figure('color','w','position',[w*(i-4),200,w,h]);i=i+1;
xdrawgraph(Xarg,yrange,method_list,'mean_reproj','Mean Reprojection For Points',...
    XLabel,'Reprojection Error (pixels)', ws);
% 
ws = ones(length(method_list), 1);
ws(1) = 0.5;
ws(3) = 0.5;
figure('color','w','position',[w*(i-4),200,w,h]);i=i+1;
xdrawgraph(Xarg,yrange,method_list,'mean_reproj_lines','Mean Reprojection For Lines',...
    XLabel,'Reprojection Error (pixels)',ws);
figure('color','w','position',[w*(i-4),200,w,h]);i=i+1;

ws = ones(length(method_list), 1);
ws(1) = 0.3;
ws(3) = 0.3;
xdrawgraph(Xarg,yrange,method_list,'mean_reproj_pts_lines','Mean Reprojection',...
    XLabel,'Reprojection Error (pixels)',ws);
% 
% figure('color','w','position',[w*(i-5),200,w,h]);i=i+1;
% xdrawgraph(Xarg,yrange,method_list,'mean_reproj_pts_lines','Mean Reprojection For Both',...
%     XLabel,'Reprojection Error (pixels)',2);

yrange= [0 0.08];
ws = ones(length(method_list), 1);
figure('color','w','position',[w*(i-4),200,w,h]);i=i+1;
xdrawgraph(Xarg,yrange,method_list,'avg_t','Time',...
    XLabel,'Average Runtime (sec)',ws);
fprintf('end');
end