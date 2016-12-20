function plot_lines_noise_5_per_row
close all;
addpath('util\');
load('lines_noise.mat');
figure('position', [100 100 1000 250]);
set(gcf,'Units','normal');

subplot(1,5,1);
XLabel = 'Noise std.dev., pix.';
Xarg = nls;
ws = ones(length(method_list), 1);
yrange = [0 1.5];
[mnames, p] = xdraw_main(Xarg,yrange,method_list,'med_r','Median Rotation',XLabel,'Rotation Error (degrees)', ws);
correct_margin();

subplot(1,5,2);
yrange = [0 7];
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'mean_r','Mean Rotation',...
    XLabel,'Rotation Error (degrees)',ws);
correct_margin();

subplot(1,5,3);
yrange= [0 2];
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'med_t','Median Translation',...
    XLabel,'Translation Error (%)',ws);
correct_margin();

subplot(1, 5, 4);
yrange= [0 5];
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'mean_t','Mean Translation',...
    XLabel,'Translation Error (%)',ws);
yrange= [0 5];
correct_margin();

subplot(1, 5, 5);
ws = ones(length(method_list), 1);
yrange = [0 10];
xdraw_main(Xarg,yrange,method_list,'mean_reproj_pts_lines','Mean Reprojection',...
    XLabel,'Reprojection Error (pixels)',ws);
correct_margin();

hL = legend(mnames, 'Orientation','horizontal', 'Position', [0.2 0.9 0.6 0.1]);


end

