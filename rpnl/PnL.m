function [rot_cw, pos_cw] = PnL(xs,xe,Vw,Pw)
%This function follows the framework of our RPnL algorithm: 
% Robust and Efficient Pose Estimation from Line Correspondences, ACCV2012.
%input: xs(:, i) = the start point of the ith image line [startpointx, startpointy, 1];
%       xe(:, i) = the end point of the ith image line   [endpointx,   endpointy, 1];
%       Vw(:, i) = the direction of ith line in the world frame
%       Pw(:, i) = a point of ith line in the world frame
%output: rot_cw = the orientation of camera in rotation matrix parametrization
%                 (V_w = rot_cw * V_c)
%        pos_cw = the position of camera in global frame;
%                 (P_w = rot_cw * P_c + pos_cw;

n = size(xs,2);
if(n~=size(xe,2) || n~=size(Vw,2) || n~=size(Pw,2)) 
    error('Input data xs, xe, Vw and Pw are inconsistent'); 
end

if(n<4)
    error('The input data is too few to determine a unique camera pose');
end

ConditionErrThreshold = 1e-3;
cosAngleThreshold     = [1.1, 0.9659, 0.8660];
optimumrot_cw = [];
optimumpos_cw = [];
lineLenVec           = zeros(n,1);
for HowToChooseFixedTwoLines = 1:3
    if HowToChooseFixedTwoLines ==1
        for i =  1:n
            lineLen = norm(xs(:,i)-xe(:,i));
            lineLenVec(i) =  lineLen;
        end
        %choose the line with longest length in the image plane;
        [longestLineLen, LineID] = max(lineLenVec);
        temp = xs(:,1);  xs(:,1) = xs(:, LineID); xs(:,LineID) = temp;
        temp = xe(:,1);  xe(:,1) = xe(:, LineID); xe(:,LineID) = temp;
        temp = Vw(:,1);  Vw(:,1) = Vw(:, LineID); Vw(:,LineID) = temp;
        temp = Pw(:,1);  Pw(:,1) = Pw(:, LineID); Pw(:,LineID) = temp;
        lineLenVec(LineID) = lineLenVec(1);
        lineLenVec(1) = 0;
        l1 = xs(:,1) - xe(:,1); l1 = l1/norm(l1);
    end
    %first line is fixed. Find the second line
    for i=2:n
        [longestLineLen, LineID] = max(lineLenVec);%the current lonest line
        l2 = xs(:,LineID) - xe(:,LineID);
        l2 = l2/norm(l2);
        lineLenVec(LineID) = 0;
        cosAngle = abs(l1'*l2);
        if cosAngle <  cosAngleThreshold(HowToChooseFixedTwoLines) % 0<angle<180, 15<angle<165,or 30<angle<150
            break;
        end 
    end
    temp = xs(:,2);  xs(:,2) = xs(:, LineID); xs(:,LineID) = temp;
    temp = xe(:,2);  xe(:,2) = xe(:, LineID); xe(:,LineID) = temp;
    temp = Vw(:,2);  Vw(:,2) = Vw(:, LineID); Vw(:,LineID) = temp;
    temp = Pw(:,2);  Pw(:,2) = Pw(:, LineID); Pw(:,LineID) = temp;
    lineLenVec(LineID) = lineLenVec(2);
    lineLenVec(2) = 0;
    
    % The rotation matrix R_wc is decomposed in way which is slightly different from the description in the paper,
    % but the framework is the same.
    % R_wc = (Rot') * R * Rot =  (Rot') * (Ry(theta) * Rz(phi) * Rx(psi)) * Rot
    
    nc1 = xcross(xs(:,1),xe(:,1));  nc1 = nc1/norm(nc1);
    Vw1 = Vw(:,1);                 Vw1 = Vw1/norm(Vw1);
    
    Xm = xcross(nc1,Vw1); Xm = Xm/norm(Xm); %the X axis of Model frame
    Ym = nc1; %the Y axis of Model frame
    Zm = xcross(Xm,Ym);  Zm= Zm/norm(Zm);%the Z axis of Model frame;
    
    Rot = [Xm, Ym, Zm]'; % Rot * [Xm, Ym, Zm] = I.
    %rotate all the vector by Rot.
    nc_bar = zeros(3,n); % nc_bar(:,i) = Rot * nc(:,i)
    Vw_bar = zeros(3,n); % Vw_bar(:,i) = Rot * Vw(:,i)
    n_c    = zeros(3,n);
    for i=1:n
        nc = xcross(xs(:,i),xe(:,i));
        nc = nc/norm(nc);
        n_c(:,i)    = nc;
        nc_bar(:,i) = Rot * nc;
        Vw_bar(:,i) = Rot * Vw(:,i);
    end
    
    
    % Determine the angle psi, it is the angle between z axis and Vw_bar(:,1).
    % The rotation matrix Rx(psi) rotates Vw_bar(:,1) to z axis
    cospsi= Vw_bar(3,1); %the angle between z axis and Vw_bar(:,1); cospsi=[0,0,1] * Vw_bar(:,1);.
    sinpsi= sqrt(1 - cospsi*cospsi);
    Rx= [1 0 0; 0 cospsi -sinpsi; 0 sinpsi cospsi];
    Zaxis = Rx * Vw_bar(:,1); % should be the Z axis, i.e. [0, 0, 1]';
    if 1 - abs(Zaxis(3)) > 1e-5
        Rx = Rx';
    end
    
    %estimate the rotation angle phi by least square residual.
    %i.e the rotation matrix Rz(phi)
    Vm2 = Rx * Vw_bar(:,2);
    A2= Vm2(1);      B2= Vm2(2);      C2= Vm2(3);
    x2= nc_bar(1,2); y2= nc_bar(2,2); z2= nc_bar(3,2);
    
    coef   = zeros(9,1); %coefficients of equation (7)
    polyDF = zeros(16,1); %dF = ployDF(1) * t^15 + ployDF(2) * t^14 + ... + ployDF(15) * t + ployDF(16);
    % construct the  polynomial F'
    for i=3:n %The first two lines are included in every triplet.
        Vm3 = Rx*Vw_bar(:,i);
        A3= Vm3(1);      B3= Vm3(2);      C3= Vm3(3);
        x3= nc_bar(1,i); y3= nc_bar(2,i); z3= nc_bar(3,i);
        u11 = -z2*A2*y3*B3 + y2*B2*z3*A3;
        u12 = -y2*A2*z3*B3 + z2*B2*y3*A3;
        u13 = -y2*B2*z3*B3 + z2*B2*y3*B3 + y2*A2*z3*A3 - z2*A2*y3*A3;
        u14 = -y2*B2*x3*C3 + x2*C2*y3*B3;
        u15 =  x2*C2*y3*A3 - y2*A2*x3*C3;
        u21 = -x2*A2*y3*B3 + y2*B2*x3*A3;
        u22 = -y2*A2*x3*B3 + x2*B2*y3*A3;
        u23 =  x2*B2*y3*B3 - y2*B2*x3*B3 - x2*A2*y3*A3 + y2*A2*x3*A3;
        u24 =  y2*B2*z3*C3 - z2*C2*y3*B3;
        u25 =  y2*A2*z3*C3 - z2*C2*y3*A3;
        u31 = -x2*A2*z3*A3 + z2*A2*x3*A3;
        u32 = -x2*B2*z3*B3 + z2*B2*x3*B3;
        u33 =  x2*A2*z3*B3 - z2*A2*x3*B3 + x2*B2*z3*A3 - z2*B2*x3*A3;
        u34 =  z2*A2*z3*C3 + x2*A2*x3*C3 - z2*C2*z3*A3 - x2*C2*x3*A3;
        u35 = -z2*B2*z3*C3 - x2*B2*x3*C3 + z2*C2*z3*B3 + x2*C2*x3*B3;
        u36 = -x2*C2*z3*C3 + z2*C2*x3*C3;
        
        a4 =   u11*u11 + u12*u12 - u13*u13 - 2*u11*u12 +   u21*u21 + u22*u22 - u23*u23...
            -2*u21*u22 - u31*u31 - u32*u32 +   u33*u33 + 2*u31*u32;
        a3 =2*(u11*u14 - u13*u15 - u12*u14 +   u21*u24 -   u23*u25...
            - u22*u24 - u31*u34 + u33*u35 +   u32*u34);
        a2 =-2*u12*u12 + u13*u13 + u14*u14 -   u15*u15 + 2*u11*u12 - 2*u22*u22 + u23*u23...
            + u24*u24 - u25*u25 +2*u21*u22+ 2*u32*u32 -   u33*u33...
            - u34*u34 + u35*u35 -2*u31*u32- 2*u31*u36 + 2*u32*u36;
        a1 =2*(u12*u14 + u13*u15 +  u22*u24 +  u23*u25 -   u32*u34 - u33*u35 - u34*u36);
        a0 =   u12*u12 + u15*u15+   u22*u22 +  u25*u25 -   u32*u32 - u35*u35 - u36*u36 - 2*u32*u36;
        b3 =2*(u11*u13 - u12*u13 +  u21*u23 -  u22*u23 -   u31*u33 + u32*u33);
        b2 =2*(u11*u15 - u12*u15 +  u13*u14 +  u21*u25 -   u22*u25 + u23*u24 - u31*u35 + u32*u35 - u33*u34);
        b1 =2*(u12*u13 + u14*u15 +  u22*u23 +  u24*u25 -   u32*u33 - u34*u35 - u33*u36);
        b0 =2*(u12*u15 + u22*u25 -  u32*u35 -  u35*u36);
        
        d0 =    a0*a0 -   b0*b0;
        d1 = 2*(a0*a1 -   b0*b1);
        d2 =    a1*a1 + 2*a0*a2 +   b0*b0 - b1*b1 - 2*b0*b2;
        d3 = 2*(a0*a3 +   a1*a2 +   b0*b1 - b1*b2 -   b0*b3);
        d4 =    a2*a2 + 2*a0*a4 + 2*a1*a3 + b1*b1 + 2*b0*b2 - b2*b2 - 2*b1*b3;
        d5 = 2*(a1*a4 +   a2*a3 +   b1*b2 + b0*b3 -   b2*b3);
        d6 =    a3*a3 + 2*a2*a4 +   b2*b2 - b3*b3 + 2*b1*b3;
        d7 = 2*(a3*a4 +   b2*b3);
        d8 =    a4*a4 +   b3*b3;
        
        coef = coef + [a4, a3, a2, a1, a0, b3, b2, b1, b0]';
        
        polyDF(1) = polyDF(1) +                                8*d8*d8;
        polyDF(2) = polyDF(2) + 15* d7*d8;
        polyDF(3) = polyDF(3) + 14* d6*d8 +                    7*d7*d7;
        polyDF(4) = polyDF(4) + 13*(d5*d8 +  d6*d7);
        polyDF(5) = polyDF(5) + 12*(d4*d8 +  d5*d7)+           6*d6*d6;
        polyDF(6) = polyDF(6) + 11*(d3*d8 +  d4*d7 +  d5*d6);
        polyDF(7) = polyDF(7) + 10*(d2*d8 +  d3*d7 +  d4*d6) + 5*d5*d5;
        polyDF(8) = polyDF(8) + 9 *(d1*d8 +  d2*d7 +  d3*d6  +   d4*d5);
        polyDF(9) = polyDF(9) + 8 *(d1*d7 +  d2*d6 +  d3*d5) + 4*d4*d4 + 8*d0*d8;
        polyDF(10)= polyDF(10)+ 7 *(d1*d6 +  d2*d5 +  d3*d4) +           7*d0*d7;
        polyDF(11)= polyDF(11)+ 6 *(d1*d5 +  d2*d4)+           3*d3*d3 + 6*d0*d6;
        polyDF(12)= polyDF(12)+ 5 *(d1*d4 +  d2*d3)+                     5*d0*d5;
        polyDF(13)= polyDF(13)+ 4 * d1*d3 +                    2*d2*d2 + 4*d0*d4;
        polyDF(14)= polyDF(14)+ 3 * d1*d2 +                              3*d0*d3;
        polyDF(15)= polyDF(15)+                                  d1*d1 + 2*d0*d2;
        polyDF(16)= polyDF(16)+                                            d0*d1;
    end
    
    %solve polyDF
    rs= roots(polyDF);
    % retriving the local minima of the cost function.
    maxreal= max(abs(real(rs)));
    rs(abs(imag(rs))/maxreal > 0.001)= [];
    minRoots = real(rs);
    
    poly    = (15:-1:1).*polyDF(1:15)';
    PolyVal = polyval(poly, minRoots);
    minRoots(PolyVal <= 0)= [];
    
    
    if isempty(minRoots)
        disp('no solution');
        return;
    end
    
    numOfRoots = length(minRoots);
    %for each minimum, we try to find a solution of the camera pose, then
    %choose the one with the least reprojection residual as the optimum of the solution.
    minimalReprojectionError = 100;
    % In general, there are two solutions which yields small re-projection error
    % or condition error:"n_c * R_wc * V_w=0". One of the solution transforms the
    % world scene behind the camera center, the other solution transforms the world
    % scene in front of camera center. While only the latter one is correct.
    % This can easily be checked by verifying their Z coordinates in the camera frame.
    % P_c(Z) must be larger than 0 if it's in front of the camera.
      
    for rootId = 1 : numOfRoots
        cosphi= minRoots(rootId);
        sign1 = sign(coef(1) * cosphi^4 + coef(2) * cosphi^3 + coef(3) * cosphi^2 + coef(4) * cosphi + coef(5));
        sign2 = sign(coef(6) * cosphi^3 + coef(7) * cosphi^2 + coef(8) * cosphi   + coef(9));
        sinphi= -sign1*sign2*sqrt(abs(1-cosphi*cosphi));
        Rz    = [cosphi -sinphi 0; sinphi cosphi 0; 0 0 1];
        %now, according to Sec4.3, we estimate the rotation angle theta
        %and the translation vector at a time.
        RzRxRot = Rz*Rx*Rot;
    
        %According to the fact that n_i^C should be orthogonal to Pi^c and Vi^c, we 
        %have: scalarproduct(Vi^c, ni^c) = 0  and scalarproduct(Pi^c, ni^c) = 0.
        %where Vi^c = Rwc * Vi^w,  Pi^c = Rwc *(Pi^w - pos_cw) = Rwc * Pi^w - pos;
        %Using the above two constraints to construct linear equation system Mat about 
        %[costheta, sintheta, tx, ty, tz, 1].
        Mat = zeros(2*n-1, 6);
        for i = 1:n
            nxi = nc_bar(1,i);  nyi = nc_bar(2,i);  nzi = nc_bar(3,i);
            Vm = RzRxRot * Vw(:,i);
            Vxi = Vm(1);        Vyi = Vm(2);        Vzi = Vm(3);
            Pm = RzRxRot * Pw(:,i);
            Pxi = Pm(1);        Pyi = Pm(2);        Pzi = Pm(3);
            % apply the constraint scalarproduct(Vi^c, ni^c) = 0
            if i>1 %if i=1, then scalarproduct(Vi^c, ni^c) always be 0
                Mat(2*i-2, 1) = nxi * Vxi + nzi * Vzi;
                Mat(2*i-2, 2) = nxi * Vzi - nzi * Vxi;
                Mat(2*i-2, 6) = nyi * Vyi;
            end
            % apply the constraint scalarproduct(Pi^c, ni^c) = 0
            Mat(2*i-1, 1) = nxi * Pxi + nzi * Pzi;
            Mat(2*i-1, 2) = nxi * Pzi - nzi * Pxi;
            Mat(2*i-1, 3) = -nxi;
            Mat(2*i-1, 4) = -nyi;
            Mat(2*i-1, 5) = -nzi;
            Mat(2*i-1, 6) = nyi * Pyi;            
        end
        %solve the linear system Mat * [costheta, sintheta, tx, ty, tz, 1]' = 0  using SVD,
        [UMat, SMat, VMat] = svd(Mat);
        vec = VMat(:,6);% the last column of Vmat;
        vec = vec/vec(6); %the condition that the last element of vec should be 1.
        normalizeTheta = 1/sqrt(vec(1)*vec(1)+vec(2)*vec(2)); %the condition costheta^2+sintheta^2 = 1;
        costheta = vec(1)*normalizeTheta;
        sintheta = vec(2)*normalizeTheta;
        Ry = [costheta, 0, sintheta; 0, 1, 0; -sintheta, 0, costheta];
                
        %now, we get the rotation matrix rot_wc and translation pos_wc
        rot_wc = (Rot') * (Ry * Rz * Rx) * Rot;
        pos_wc = -Rot'*vec(3:5);     
        
        %now normalize the camera pose by 3D alignment. We first translate the points
        %on line in the world frame Pw to points in the camera frame Pc. Then we project
        %Pc onto the line interpretation plane as Pc_new. So we could call the point
        %alignment algorithm to normalize the camera by aligning Pc_new and Pw.
        %In order to improve the accuracy of the aligment step, we choose two points for each
        %lines. The first point is Pwi, the second point is  the closest point on line i to camera center. 
        %(Pw2i = Pwi - (Pwi'*Vwi)*Vwi.) 
        Pw2    = zeros(3,n);
        Pc_new = zeros(3,n);
        Pc2_new = zeros(3,n);
        for i=1:n
            nci = n_c(:,i);
            Pwi = Pw(:,i);
            Vwi = Vw(:,i);
            %first point on line i
            Pci     = rot_wc*Pwi + pos_wc;
            Pc_new(:,i) = Pci - (Pci'*nci)*nci;
            %second point is the closest point on line i to camera center.
            Pw2i    = Pwi - (Pwi'*Vwi)*Vwi;
            Pw2(:,i)= Pw2i;
            Pc2i    = rot_wc*Pw2i + pos_wc;
            Pc2_new(:,i) = Pc2i - (Pc2i'*nci)*nci;
        end

        [rot_wc,pos_wc] = calcampose([Pc_new Pc2_new],[Pw Pw2]);
        pos_cw = -rot_wc'*pos_wc;
        
        
        %check the condition n_c^T * rot_wc * V_w = 0;
        conditionErr = 0;
        for i=1:n
            conditionErr = conditionErr + (n_c(:,i)' * rot_wc * Vw(:,i))^2;
        end
        
        if conditionErr/n < ConditionErrThreshold || HowToChooseFixedTwoLines ==3
            %check whether the world scene is in front of the camera.
            numLineInFrontofCamera = 0;
            if HowToChooseFixedTwoLines <3
                for i=1:n
                    P_c = rot_wc*(Pw(:,i) - pos_cw);
                    if P_c(3) > 0
                        numLineInFrontofCamera = numLineInFrontofCamera+1;
                    end
                end
            else 
                numLineInFrontofCamera = n;
            end
            
            if numLineInFrontofCamera>0.5*n
                %most of the lines are in front of camera, then check the reprojection error.
                reprojectionError = 0;
                for i=1:n
                    nc = rot_wc * xcross(Pw(:,i) - pos_cw, Vw(:,i));%line projection function
                    h1 = nc' * xs(:,i);
                    h2 = nc' * xe(:,i);
                    lineLen = norm(xs(:,i)-xe(:,i))/3;
                    reprojectionError = reprojectionError +  lineLen * (h1*h1 + h1*h2 + h2*h2) / (nc(1)*nc(1)+nc(2)*nc(2));
                end
                if reprojectionError < minimalReprojectionError
                    optimumrot_cw = rot_wc';
                    optimumpos_cw = pos_cw;
                    minimalReprojectionError = reprojectionError;
                end
            end
        end
    end    
     
    if  size(optimumrot_cw,1)>0 
        break;
    end
end

rot_cw = optimumrot_cw;
pos_cw = optimumpos_cw;

return

function c = xcross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
 
 return