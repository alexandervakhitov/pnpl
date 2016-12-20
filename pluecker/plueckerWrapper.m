function [R,t] = plueckerWrapper(xs, xe, Xs, Xe)
    lineNum = size(xs, 2);
    line_3D_end_pts = ones(4, 2*lineNum);
    line_2D_end_pts = ones(3, 2*lineNum);
    for i = 1:lineNum
        line_3D_end_pts(1:3, 2*i-1) = Xs(:, i);
        line_3D_end_pts(1:3, 2*i) = Xe(:, i);
        line_2D_end_pts (1:2, 2*i-1) = xs(:, i);
        line_2D_end_pts (1:2, 2*i) = xe(:, i);
    end
    w = ones(lineNum, 1);
    line_2D_end_pts(2,:) = -line_2D_end_pts(2,:);
    [R,t] = linePoseEstim( line_3D_end_pts, line_2D_end_pts, w );
    t = R*t;
    R = diag([-1 1 -1])*R;            
    t = diag([-1 1 -1])*t;
end