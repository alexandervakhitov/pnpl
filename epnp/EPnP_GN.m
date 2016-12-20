function [R,t] = EPnP_GN(XX,xx)

[R,t]= epnp_orig.efficient_pnp_gauss(XX.',xx.',diag([1 1 1]));

return