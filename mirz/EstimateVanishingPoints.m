function [q_sols, residuals, info] = EstimateVanishingPoints(I_moment_n, G_line, method, opts)

% Estimate Vanishing Points for a Calibrated Camera
% $Id: EstimateVanishingPoints.m 1262 2012-01-20 01:55:00Z faraz $

% measurements
n = I_moment_n;
e = G_line;
N = size(n, 2);

% method
%relaxed = strcmp(config,'relaxed');
%assert(relaxed || strcmp(config,'nonrelaxed'), 'The problem can be solved
%either in relaxed form or nonrelaxed form.')
method = lower(method);

if nargin < 4
    % generate default options
    opts = EstimateVanishingPointsDefaultOpts(method);
end

conditioning = opts.conditioning;
outputFormat = opts.outputFormat;

kActionFunctionIndex = opts.kActionFunctionIndex;
kNormalizationTol = opts.kNormalizationTol;
kPolyResidualTol = opts.kPolyResidualTol;
kSolutionInClassTol = opts.kSolutionInClassTol;
kMulMatrixCondToAbort = opts.kMulMatrixCondToAbort;
Ne = opts.Ne;
NeMax = opts.NeMax;


% conditioning
if conditioning
    q_trans = [rand(3,1)*2-1; rand(1)];
    q_trans = q_trans/norm(q_trans);
else
    q_trans = [0; 0; 0; 1];
end
e_trans = quat2rot(q_trans)*e;



% load the polynomial and normal set exponents and construct the polynomials
switch method
    case 'minimal'
        NormalSetExp = GenNormalSetMinimal();
        p = GetPolysMinimal(n,e_trans); % uses only the first 3 pairs
    case 'relaxed'
        NormalSetExp = GenNormalSetRelaxed();
        p = GetPolysRelaxed(n,e_trans,N);
    case 'nonrelaxed'
        NormalSetExp = GenNormalSet();
        p = GetPolys(n,e_trans,N);
    otherwise
        error('Uknown method requested for estimating vanishing points.');
end



% normalization
for k = 1:length(p)
    p(k).coef = p(k).coef/max(abs(p(k).coef));
end



% action function
act_fun = zeros(size(NormalSetExp,1),1);
act_fun(kActionFunctionIndex) = 1;

% compute the multiplication matrix
% increase Ne until desirable cond number is acheived
tic;
info.mulMatrixCompCond = inf;
while (info.mulMatrixCompCond > kMulMatrixCondToAbort && Ne <= NeMax)
    [m_f, mulMatrixInfo] = GetMulMatrixCash(NormalSetExp, p, act_fun, Ne, kMulMatrixCondToAbort);
    Ne = Ne+1;
end
info.mulMatrixCompTime = toc;
info.mulMatrixCompCond = mulMatrixInfo.cond;

% if the multiplication's matrix computation involved inverting a matrix
% with condition number higher than the given threshold, no reliable
% solution can beextracted. Abort.
if info.mulMatrixCompCond > kMulMatrixCondToAbort
    q_sols = []; 
    residuals = []; 
    return
end



% eigenvalue decomposition of the multiplication matrix
[V ~] = eig(m_f');


% the last element corresponding to a single variable
switch method
    case {'relaxed','minimal'}
        lastVariableIndex = 4;
    case 'nonrelaxed'
        lastVariableIndex = 5;
end

% remove unstable eigenvalues that have first element too close to zero
% also remove complex eigenvalues
stableRealEig = (abs(V(1,:)) > kNormalizationTol) & ...
                ~any(imag(V(lastVariableIndex:-1:1,:)));
% if no stable and real eigenvector, abort            
if sum(stableRealEig) == 0
    q_sols = []; 
    residuals = []; 
    return
end

% normalize the eigenvalues to             
polySols = V(lastVariableIndex:-1:2,stableRealEig)./ ...
           (ones(lastVariableIndex-1,1)*V(1,stableRealEig));
       
       

% check the polynomial residuals
polyResisual = zeros(1,size(polySols,2));
for k = 1:size(polyResisual,2)
    polyResisual(1,k) = norm(PolyEval(p, polySols(:,k)));
end

[polyResisual ind] = sort(polyResisual);
polySols = polySols(:,ind);
validMinimums = (polyResisual < kPolyResidualTol);


if strcmp(method, {'relaxed','nonrelaxed'})
    % Compute Hessian and check if it is PSD (note that H is symmetric and only
    % upper triangle is computed
    hessian = PolyHessian(p);
    %ss = [];
    for k = 1:size(polySols,2)
        hessianValue = PolyEval(hessian,polySols(:,k));
        % make it symmetric
        hessianValue = hessianValue + triu(hessianValue,1)';
        [~,isNotPSD] = chol(hessianValue);
        %ss = [ss eig(hessianValue)];
        if isNotPSD == false
            validMinimums(k) = 1; % should we touch this? probably not.
        else
            validMinimums(k) = 0;
        end
    end
end


% compute quaternion, and sort them according to polyResiduals
% validResiduals = polyResisual(validMinimums);
q_sols = quat_mul_batch(quat_inv(q_trans),cgr2quat(polySols(1:3,validMinimums)));


% if enabled, replace class solutions with the best one
if outputFormat ~= 0
    q_sols_proc = [];
    processedSols = false(1,size(q_sols,2));
    for i = 1:length(processedSols)
        % if the solution is already processed, skip it
        if (processedSols(i) == false)
            % otherwise add it (or its expansion) to the set;
            if outputFormat == 1
                q_sols_proc = [q_sols_proc createClassSolution(q_sols(:,i))];
            elseif outputFormat == 2
                q_sols_proc = [q_sols_proc q_sols(:,i)];
            else
                error('output format not recognized!')
            end
            
            for j = i+1:length(processedSols)
                if norm(quat_mul(q_sols(:,i), quat_inv(q_sols(:,j))),inf) > kSolutionInClassTol
                    processedSols(j) = true;
                end
            end
        end
    end
    q_sols = q_sols_proc;
end


residuals = zeros(1,size(q_sols,2));
for k = 1:size(q_sols,2)
    for j = 1:N
        residuals(k) = residuals(k) + (e(:,j)'*quat2rot(q_sols(:,k))*n(:,j))^2;
    end
end


[residuals, ind] = sort(residuals);
q_sols = q_sols(:,ind);


return
