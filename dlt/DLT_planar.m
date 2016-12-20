function [R,t] = DLT_planar(Xw,xx,xs, xe, Xs, Xe, lnc)
    [R, t] = DLT_main(Xw(1:2,:),xx,xs, xe, Xs(1:2,:), Xe(1:2,:), 1);    
end