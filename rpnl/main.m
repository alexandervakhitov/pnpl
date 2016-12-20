%test the pose estimation from line correspondences problem
clear all;
clc;

numOfLines = 5; %please set the number of lines (numOfLines>=4)
noiseStd   = 5/1024; %please set the Std of image noise. (noiseStd = 0 or 10 pixels)

%first creat simulation data; 
startPoint =  rand(numOfLines,2) - 0.5;%creat the end points of line in the image plane
endPoint   =  rand(numOfLines,2) - 0.5;

figure;
hold on;
for i=1:numOfLines
plot(startPoint(i,:), endPoint(i,:))
end


startPointDepth = 2*rand(numOfLines,1)+1.1;%creat the depth of line end points to generate 3D line end points
endPointDepth   = 2*rand(numOfLines,1)+1.1;
startPoint3D = [startPoint(:,1).* startPointDepth, startPoint(:,2).* startPointDepth,startPointDepth];
endPoint3D = [endPoint(:,1).* endPointDepth, endPoint(:,2).* endPointDepth, endPointDepth];

%Rotation matrix from camera to world frame R_cw
theta_x = 2*rand(1)-1;
theta_y = 2*rand(1)-1;
theta_z = 2*rand(1)-1;
r_x = [1,0,0; 0, cos(theta_x), - sin(theta_x); 0, sin(theta_x), cos(theta_x) ];
r_y = [cos(theta_y), 0, sin(theta_y); 0, 1, 0; - sin(theta_y), 0, cos(theta_y)];
r_z = [cos(theta_z), - sin(theta_z), 0; sin(theta_z), cos(theta_z), 0; 0, 0, 1];
R_cw = r_z * r_y * r_x
% translation from camera to world frame T_cw
T_cw = 2*rand(3,1)

%generate the 3D line direction V and point P in world frame
startPoint3D_w = zeros(numOfLines,3);
endPoint3D_w   = zeros(numOfLines,3);
V = zeros(numOfLines,3);
for i=1:numOfLines 
    startPoint3D_w(i,:) = (R_cw * startPoint3D(i,:)' + T_cw)'; %P_w = R_cw * P_c + T_cw i.e. P_c = R_wc (P_w - T_cw);
    endPoint3D_w(i,:)   = (R_cw * endPoint3D(i,:)' + T_cw)';
    dir = R_cw*(endPoint3D(i,:) - startPoint3D(i,:))';
    V(i,:) = dir/norm(dir);
end
P = startPoint3D_w;
%add some noise to the endpoints 
startPointNoise =  noiseStd*(rand(numOfLines,2) - 0.5);
endPointNoise   =  noiseStd*(rand(numOfLines,2) - 0.5);
startPoint      =  startPoint + startPointNoise;
endPoint        =  endPoint   + endPointNoise;

xs = zeros(3,numOfLines);
xe = zeros(3,numOfLines);
for i= 1:numOfLines
    xs(1:3, i) = [startPoint(i,:), 1]';
    xe(1:3, i) = [endPoint(i,:), 1]';
end

for i=1:numOfLines
plot(startPoint(i,:), endPoint(i,:), 'r')
end
hold off;
 
%call the PnL() function to estimate camera pose.
[RPnLR_cw, RPnLT_cw] = PnL(xs,xe,V',P')
%call the R_and_T() function to refine the estimated camera pose.
[RPnLppR_cw, RPnLppT_cw] =  R_and_T(xs,xe,startPoint3D_w',endPoint3D_w',RPnLR_cw, RPnLT_cw)
