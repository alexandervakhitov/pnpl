function H = PolyHessian(p)

% $Id: PolyHessian.m 1264 2012-01-20 02:03:57Z faraz $

varNum = size(p(1).exp,1);

H(varNum,varNum) = struct('exp',[],'coef',[],'m',[]);

for i = 1:varNum
    for j = i:varNum
        
        H(i,j) = derivative(p(i),j);
        
    end
end



function d = derivative(d,j)

coef = d.exp(j,:).*d.coef;
d.exp = d.exp(:,coef~=0);
d.exp(j,:) = d.exp(j,:) - 1;
d.coef = coef(coef~=0);
d.m = length(d.coef);


