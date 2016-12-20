function [R, t] = PNLWrapper(xs, xe, Xs, Xe)
    P = Xs;
    V = Xe-Xs;
    for i = 1:size(V, 2)
       V(:, i) = V(:, i) / norm(V(:, i)); 
    end
    xs = [xs; ones(1, size(xs, 2))];
    xe = [xe; ones(1, size(xe, 2))];
    [RPnLR_cw, RPnLT_cw] = PnL(xs,xe,V,P);
    if (size(RPnLR_cw, 1) < 3)
        R = [];
        t = [];
        return;
    end
    [RPnLppR_cw, RPnLppT_cw] =  R_and_T(xs,xe,Xs,Xe,RPnLR_cw, RPnLT_cw);
    R = RPnLR_cw';
    t = -R*RPnLT_cw;
end