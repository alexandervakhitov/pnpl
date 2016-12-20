function [R,t] = EPnP_GNN(XX,xx)

[R,t]= epnp_orig.efficient_pnp_gauss_new(XX.',xx.',diag([1 1 1]));

return