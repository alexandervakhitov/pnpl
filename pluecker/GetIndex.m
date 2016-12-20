function ind = GetIndex(exp)

% get the exponent vector and return the index, returns -1 if exp does not
% exist
% $Id: GetIndex.m 1264 2012-01-20 02:03:57Z faraz $

if any( exp < 0 ) 
    ind = -1;
    return
end

N = length(exp);
d = sum(exp);

% greater than or equal index
ge_ind = 0;
% exponents that are already accounted for
exp_pre = 0;
for k = N:-1:1
    if exp(k) == 0
        continue;
    end
        
    for digit = 0:exp(k)-1
        %nchoosek((k-1)+(d-digit-exp_pre)-1, (k-1)-1)
        if (k-1)-1 >= 0
            ge_ind = ge_ind + nchoosek((k-1)+(d-digit-exp_pre)-1, (k-1)-1);
        end
    end
    exp_pre = exp_pre + exp(k);
end

ind = nchoosek(N+d, N) - ge_ind ;

return
