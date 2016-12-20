function [Rot_cw, Pos_cw] = R_and_T(xs,xe, P1w, P2w, initRot_cw, initPos_cw, maxIterNum, TerminateTh)
%This function follows the iterative algorithm (R_and_T) proposed by Kumar and Hanson: 
% Robust methods for estimating pose and a sensitivity analysis.
%input: xs(:, i)  = the start point of the ith image line [startpointx, startpointy, 1];
%       xe(:, i)  = the end point of the ith image line   [endpointx,   endpointy, 1];
%       P1w(:, i) = one point of ith line in the world frame, not nessary to be the endpoints.
%       P2w(:, i) = another point of ith line in the world frame, not nessary to be the endpoints.
%       initRot_cw: the iterative algorithm need initial value of the roation matrix
%       initPos_cw: the iterative algorithm need initial value of the translation vector 
%       maxIterNum: the maximum number of iterations allowed by the  iterative algorithm
%       TerminateTh:the threshold to accept the estimated results and termiate the algorithm.
%output: Rot_cw = the orientation of camera in rotation matrix parametrization
%                 (V_w = Rot_cw * V_c)
%        Pos_cw = the position of camera in global frame; 
%                 (P_w = Rot_cw * P_c + pos_cw;
%                  P_c = Rot_wc * P_w + pos_wc;
%                  Pos_wc = -Rot_cw' * Pos_cw and Pos_cw = -Rot_wc' * Pos_wc)

if nargin < 6
    error('Please do not forget the initial value.');
end

if nargin < 8
    % set the default paprameters
    maxIterNum  = 20;
    TerminateTh = 1e-5; 
end

n = size(xs,2);
if(n~=size(xe,2) || n~=size(P1w,2) || n~=size(P2w,2)) 
    error('Input data xs, xe, P1w and P2w are inconsistent'); 
end

if(n<4)
    error('The input data is too few to determine a unique camera pose');
end

%first compute the weight of each line and the normal of the interpretation plane passing through to camera center and the line 
w  = zeros(n,1);
nc = zeros(3,n);
for i=1:n
    w(i) = 1/norm(xs(:,i)-xe(:,i)); % the weight of a line is the inverse of its image length
    temp = cross(xs(:,i),xe(:,i));  
    nc(:,i) = temp/norm(temp);   
end

Rot_wc = initRot_cw';
Pos_wc = -initRot_cw' * initPos_cw;
for iter = 1: maxIterNum
    %construct the equation (31)
    A = zeros(6,7); % A*[dT, dOmiga, 1] = 0
    C = zeros(3,3);
    D = zeros(3,3);
    F = zeros(3,3);
    c_bar = zeros(3,1);
    d_bar = zeros(3,1);
    for i = 1:n
        Pi = Rot_wc * P1w(:,i);% for first point on line
        Ni = nc(:,i);
        wi = w(i);
        bi = cross(Pi, Ni);
        C  = C + wi*Ni*Ni';
        D  = D + wi*bi*bi';
        F  = F + wi*Ni*bi';
        scale = wi*(Ni'*(Pi+Pos_wc));
        c_bar = c_bar + scale*Ni;
        d_bar = d_bar + scale*bi;
        Pi = Rot_wc * P2w(:,i);% for second point on line
        Ni = nc(:,i);
        wi = w(i);
        bi = cross(Pi, Ni);
        C  = C + wi*Ni*Ni';
        D  = D + wi*bi*bi';
        F  = F + wi*Ni*bi';
        scale = wi*(Ni'*(Pi+Pos_wc));
        c_bar = c_bar + scale*Ni;
        d_bar = d_bar + scale*bi;
    end
    A(1:3, 1:3) = C;    A(1:3, 4:6) = F;    A(1:3, 7) = c_bar;
    A(4:6, 1:3) = F';   A(4:6, 4:6) = D;    A(4:6, 7) = d_bar;
    %sovle the system by using SVD;
    [UMat, SMat, VMat] = svd(A);
    vec = VMat(:,7);% the last column of Vmat;
    vec = vec/vec(7); %the condition that the last element of vec should be 1.
    
    %update the rotation and translation parameters;
    dT = vec(1:3);  dOmiga =  vec(4:6);    
    Rot_wc = [1, -dOmiga(3), dOmiga(2); dOmiga(3), 1, -dOmiga(1); -dOmiga(2), dOmiga(1), 1]*Rot_wc; % newRot_wc = ( I + [dOmiga]x ) oldRot_wc
	% !!!!!!!! may be we can compute new R using rodrigues(r+dr)
    Pos_wc = Pos_wc + dT;
	

    if norm(dT)<TerminateTh && norm(dOmiga)<0.1*TerminateTh
        break;
    end
end

Rot_cw = Rot_wc';
Pos_cw = -Rot_cw*Pos_wc;


return