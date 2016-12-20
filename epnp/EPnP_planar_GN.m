function [R,t] = EPnP_planar_GN(XX,xx)

    [R,t]= epnp_orig.efficient_pnp_planar_gauss(XX.',xx.',diag([1 1 1]));

return