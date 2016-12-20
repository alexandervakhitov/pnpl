function f = PolyEval(p,s)

% $Id: PolyEval.m 1264 2012-01-20 02:03:57Z faraz $

f = zeros(size(p));

%assert(size(s,2)==1,'solution must be a column vector')

for i = 1:size(p,1)
    for j = 1:size(p,2)
        
        if isempty(p(i,j).m) || p(i,j).m == 0
            f(i,j) = 0;
        else
            f(i,j) = prod((s*ones(1,p(i,j).m)).^p(i,j).exp)*p(i,j).coef';
        end
        
    end
end
