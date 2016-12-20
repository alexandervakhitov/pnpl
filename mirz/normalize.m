function y = normalize(x, fc, cc, kc, alpha_c)

assert(kc == 0, 'Nonlinear distortion removal is not implemented. Use CalTech Camera Calibration Toolbox')
assert(alpha_c == 0, 'Nonzero skew is not implemented. Use CalTech Camera Calibration Toolbox')

y = (x - cc)./fc;
