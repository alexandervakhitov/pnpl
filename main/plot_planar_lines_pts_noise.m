function plot_lines_pts_noise
load('planar_lines_pts_noise.mat');

close all;
addpath('util\');
h = sp_position();
set(gcf,'Units','normal');
set(gcf, 'PaperPositionMode', 'auto');

% subplot(1,5,1);
sp_format(1);
XLabel = 'Noise, std.dev. (pix.)';
Xarg = nls;
ws = ones(length(method_list), 1);
yrange = [0 1.5];
[mnames, p] = xdraw_main(Xarg,yrange,method_list,'med_r','Median Rotation',XLabel,'Rotation Error (degrees)', ws);
correct_margin();

% subplot(1,5,2);
sp_format(2);
yrange = [0 7];
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'mean_r','Mean Rotation',...
    XLabel,'Rotation Error (degrees)',ws);
correct_margin();

% subplot(1,5,3);
sp_format(3);
yrange= [0 2];
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'med_t','Median Translation',...
    XLabel,'Translation Error (%)',ws);
correct_margin();

% subplot(1, 5, 4);
sp_format(4);
yrange= [0 5];
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'mean_t','Mean Translation',...
    XLabel,'Translation Error (%)',ws);
yrange= [0 5];
correct_margin();

% subplot(1, 5, 5);
% sp_format(5);
% ws = ones(length(method_list), 1);
% yrange = [0 10];
% xdraw_main(Xarg,yrange,method_list,'mean_reproj_pts_lines','Mean Reprojection',...
%     XLabel,'Reprojection Error (pixels)',ws);
% correct_margin();

% hL = legend(mnames, 'Orientation','horizontal', 'Position', [0.2 0.9 0.6 0.1]);
hL = legend(mnames, 'Orientation','vertical', 'Position', [0.85 0.2 0.1 0.6]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures\filename.pdf','-dpdf','-r0')
end