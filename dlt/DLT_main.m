function [R,t] = DLT_main(Xw,xx,xs, xe, Xs, Xe, is_planar)
    npts = size(Xw, 2);        
    x = xx(1, :)';
    y = xx(2, :)';
    nl = size(xs, 2);
    l2d = epnp_lines.pnl_preprocess(xs, xe);
    a = l2d(1, :)';
    b = l2d(2, :)';
    c = l2d(3, :)';    
    if (is_planar == 1)
        A_pts_1 = [zeros(npts, 2), -Xw', repmat(y, 1, 2).*Xw', zeros(npts, 1), -ones(npts, 1), y];
        A_pts_2 = [Xw', zeros(npts, 2), -repmat(x, 1, 2).*Xw', ones(npts, 1), zeros(npts, 1), -x];            
        A_lns_1 = [repmat(a, 1, 2).*Xs', repmat(b, 1, 2).*Xs', repmat(c, 1, 2).*Xs', a, b, c];
        A_lns_2 = [repmat(a, 1, 2).*Xe', repmat(b, 1, 2).*Xe', repmat(c, 1, 2).*Xe', a, b, c];        
    else
        A_pts_1 = [zeros(npts, 3), -Xw', repmat(y, 1, 3).*Xw', zeros(npts, 1), -ones(npts, 1), y];
        A_pts_2 = [Xw', zeros(npts, 3), -repmat(x, 1, 3).*Xw', ones(npts, 1), zeros(npts, 1), -x];            
        A_lns_1 = [repmat(a, 1, 3).*Xs', repmat(b, 1, 3).*Xs', repmat(c, 1, 3).*Xs', a, b, c];
        A_lns_2 = [repmat(a, 1, 3).*Xe', repmat(b, 1, 3).*Xe', repmat(c, 1, 3).*Xe', a, b, c];
    end
    
    A = [A_pts_1; A_pts_2; A_lns_1; A_lns_2];
    
    global Rtt;
    global tt;
    
    [U,S,V] = svd(A);
    sol = V(:, end);
    if (is_planar == 1)
        P = [reshape(sol(1:6), 2, 3)', sol(7:9)];
        pvec = P(:);
    else        
        P = [reshape(sol(1:9), 3, 3)', sol(10:12)];
        pvec = P(:);
    end
    
    X = [Xw Xs Xe];
    Xh = [X; ones(1, size(X, 2))];
    Alph = Xh';
    
    [Cc,Xc,sc]=compute_norm_sign_scaling_factor(pvec,Alph,X');
    if (is_planar == 1)
        X = [X; zeros(1, size(X, 2))];
%         Xc = [Xc zeros(size(Xc, 1), 1)];
    end
    [R,t]=getrotT(X',Xc); 
    if (is_planar == 1)
        R(:, 3) = cross(R(:, 1), R(:, 2));        
    end
end