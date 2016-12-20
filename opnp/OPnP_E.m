function [R,t] = OPnP_E(xs, xe, Xs, Xe, lnc)
    xx = [xs xe];
    Xw = [Xs Xe];
    [R,t] = OPnP(Xw, xx);
end
