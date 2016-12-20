function [m_f, info] = GetMulMatrixCash(rcoeffs, p, act_fun, Ne, kMulMatrixCondToAbort)

% GetMulMatrixCash computes the multiplication matrix, from a given set of polynomials,
% their normal set, action (multiplication) function, and expansion order. The first time
% that this function runs, it computes the structure of Cexp, and cashes it for future use. 
% Therefore, the first invocation for a specific normal set, action function, and exp order
% may take considerably longer than consequent runs.
%
% Input arguments: 
%
% rcoeffs: A matrix encoding exponents of the monomial in the normal set. Number of rows 
% of this matrix corresponds to the number of variables, and number of columns equals
% the cardinality of normal set (i.e., number of solution). 
% 
% p: Array of structure encoding the polynomials, with three fields:
% p().exp: exponents of the monomials
% p().coef: coefficients corresponding to monomials
% p().m: number of monomials involved
%
% Ne: Expansion order
%
% kMulMatrixCondToAbort: maximum acceptable condition number of unreduced mult. matrix.
% if condition number of unreduced multiplication matrix is higher than 
% kMulMatrixCondToAbort, the multiplication matrix cannot be computed 
% reliably. The function will return m_f = [], and aborts.
%
% This function uses DataHash.m function developed and copyrighted by Jan Simon.
% Refer to DataHash-license.txt regarding licensing of DataHash.m
%
% Author: Faraz Mirzaei, University of Minnesota. 
% contact faraz -at- umn.edu
% Copyright (c) 2010-2012 The Regents of the University of Minnesota
%
% $Id: GetMulMatrixCash.m 1376 2012-07-27 05:21:36Z faraz $

% number of polynomials
Np = size(p,2);

% number of variables
Nv = size(rcoeffs,1);

% total number of monomials up to and including degree 9
N_monom = nchoosek(Nv+Ne,Nv);

% highest total degree monomial in normal set
db = max(sum(rcoeffs));
assert(Ne > max(db), 'Ne needs to be large enough to cover all xr monomials');

% total degree of polynomial (assuming all have the same degree)
dp = zeros(1,Np);
for i = 1:Np
    dp(i) = max(sum(p(i).exp));
end

% number of monomials in normal set xb (i.e., # of solutions)
lb = size(rcoeffs,2);

% number of extender polynomials applicable to each polynomial
N_ext = zeros(1,Np);
for i = 1:Np
    N_ext(i) = nchoosek(Nv+Ne-dp(i), Nv);
end

% number of nonzero elements in Cexp
Cexp_nnz = [p.m]*N_ext';

identifier = DataHash({rcoeffs,act_fun,p.exp,Ne});
cash_filename = strcat('cash/',identifier,'.mat');
% check the cash
regenerate_Cexp_structure = true;
if exist(cash_filename,'file')
    load(cash_filename);
    if size(Cexp_structure,1) == Cexp_nnz ...
            && isequal(check_act_fun,act_fun) ...
            && isequal(check_rcoeffs, rcoeffs) %#ok<NODEF>
        regenerate_Cexp_structure = false;
    end
end


% regenerate variables if cash not available
if regenerate_Cexp_structure == true

    fprintf('cash file does not exist, or its does not match the given system. \n')
    fprintf('It will be generated and cashed now. This may take a while ...\n')
   
    % The total number of extender monomials that have to be generated
    % Note that not all polynomials use all the extender monomials for
    % expansion
    N_ext_all = nchoosek(Nv+Ne-min(dp), Nv);
    extender = zeros(Nv, N_ext_all);
    % extender include neutral element:
    % i.e. extender(:,1)*p(j) is polynomial p(j) itself
    extender(:,1) = zeros(Nv,1);
    for k = 2:N_ext_all
        extender(:,k) = NextDeg2(extender(:,k-1));
    end
    

    ind_map = zeros(1,N_monom);
    
    
    % fill in ind_map for basis monomials (xb)
    for k = 1:lb
        ind = GetIndex(rcoeffs(:,k));
        assert(ind_map(ind) == 0, 'double indexing detected');
        ind_map(ind) = k;
    end
    
    % fill in ind_map for xr and find the size of xr (lr)
    lr = 0;
    for k = 1:lb
        exp = rcoeffs(:,k) + act_fun;
        ind = GetIndex(exp);
        if ind_map(ind) == 0
            lr = lr + 1;
            ind_map(ind) = lb + lr;
            %exp
        end
    end
    
    % pointers to start of xr and xe
    eptr = lb + lr;
        
    Cexp_cntr = 1;
    Cexp_structure = zeros(Cexp_nnz,4);
    
    % make this a loop for all j
    for j = 1:Np
        % for all extenders
        for k = 1:N_ext(j)
            % for all monomials of plynomial p(j)
            for l = 1:p(j).m
                % extend the monomial
                exp = extender(:,k) + p(j).exp(:,l);
                
                % get the index of the extended monomial
                ind = GetIndex(exp);
                
                %lookup monomial position in C
                mon_ptr = ind_map(ind);
                
                % if mon_ptr is not set yet, it should belong to xe:
                if mon_ptr == 0
                    
                    eptr = eptr + 1;
                    ind_map(ind) = eptr;
                    
                    % update mon_ptr to reflect the new assignment
                    mon_ptr = ind_map(ind);
                    
                end
                
                %Cexp_i(Cexp_cntr) = sum(N_ext(1:j-1))+k;
                %Cexp_j(Cexp_cntr) = mon_ptr;
                %Cexp_val(Cexp_cntr) = p(j).coef(l);
                Cexp_structure(Cexp_cntr,:) = [sum(N_ext(1:j-1))+k, ...
                    mon_ptr, ...
                    j, ...
                    l];
                Cexp_cntr = Cexp_cntr+1;
                
            end
        end
    end
    
    check_act_fun = act_fun; %#ok<NASGU>
    check_rcoeffs = rcoeffs; %#ok<NASGU>
    [~, ~, ~] = mkdir('cash');
    save(cash_filename, 'Cexp_structure', 'lr', 'eptr', 'ind_map', ...
        'extender', 'check_act_fun', 'check_rcoeffs');
    fprintf('cash file generated and saved now.\n')
    
end

Cexp_val = zeros(Cexp_nnz,1);
for k = 1:Cexp_nnz
    Cexp_val(k) = p(Cexp_structure(k,3)).coef(Cexp_structure(k,4));
end

Cexp = sparse(Cexp_structure(:,1), Cexp_structure(:,2), Cexp_val, sum(N_ext), N_monom);
info.CexpSize = [sum(N_ext), N_monom];
info.CexpPartitionInfo = [lb, lr, eptr-lb-lr];
info.CexpNNZ = Cexp_nnz;

tol = 1e-8; % tolerance for finding nullspace

CB = Cexp(:,1:lb);
CR = Cexp(:,lb+1:lb+lr);
CE = Cexp(:,lb+lr+1:eptr);



%%{
% sparse version
[SP_QtCRCB, SP_RCE, ~] = qr(CE, [CR CB]); % the lower portion of QtCRCB has the left nullspace premultiplied to the rest. Find out the rows next
%SP_nullinds = find(full(max(abs(SP_RCE),[],2)> tol),1,'last')+1 : size(Cexp,1); % rows of nullspace
nullBegin = find(abs(diag(SP_RCE)) > tol,1,'last')+1;
%info.cond = cond(full(SP_QtCRCB(nullBegin:end, 1:lr)));
%T = SP_QtCRCB(nullBegin:end, 1:lr) \ SP_QtCRCB(nullBegin:end, lr+1:end);
[C, R] = qr(SP_QtCRCB(nullBegin:end, 1:lr),SP_QtCRCB(nullBegin:end, lr+1:end),0);

%CONDEST invokes RAND. so save the state and reload it

% ipribyl CHANGE --- FROM (update from obsolete random number generator)
% saveState = rand('state');
% ipribyl CHANGE -- -TO
saveState = rng;
% ipribyl CHANGE --- END

if size(R,1) ~= size(R,2)
   info.cond = inf;
else
    info.cond = condest(R);
end    

% ipribyl CHANGE --- FROM (update from obsolete random number generator)
% rand('state',saveState);
% ipribyl CHANGE -- -TO
rng(saveState);
% ipribyl CHANGE --- END

% if condition number is too high, abort
if info.cond > kMulMatrixCondToAbort
    m_f = [];
    return
end

T=R\C;
%%}


%{
[QCE, RCE, e] = qr(full(CE)); % qr with column permutations. more stable
% toc
clear e
nullinds = find(full(max(abs(RCE),[],2)> tol),1,'last')+1 : size(Cexp,1); % rows of nullspace
clear RCE
QtCRCB = QCE(:,nullinds)'*[CR CB]; % need to do this by hand :-(
clear QCE
fprintf('size of CR: %dx%d\nrank of CR: %d \n',size(CR,1),size(CR,2),rank(QtCRCB(:, 1:lr)));
T = QtCRCB(:, 1:lr) \ QtCRCB(:, lr+1:end);
clear QtCRCB
%}

% CEperp=null(full(CE)');
% NtCRCB = CEperp'*[CR CB];
% T = NtCRCB(:, 1:lr) \ NtCRCB(:, lr+1:end);


% unreduced multiplication matrix
Mxk = zeros(lr+lb, lb);
for k = 1:lb
    exp = rcoeffs(:,k) + act_fun;
    ind = GetIndex(exp);
    assert(ind_map(ind) ~= 0,'unexpected zero in ind_map')
    Mxk(ind_map(ind), k) = 1;
end

P = [eye(lb) -T'];

m_f = full(P*Mxk);









