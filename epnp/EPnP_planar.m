function [R,t] = EPnP_planar(XX,xx)

[R,t]= efficient_pnp_planar(XX.',xx.',diag([1 1 1]));
%     [R,t]= efficient_pnp_planar_gauss(XX.',xx.',diag([1 1 1]));

return