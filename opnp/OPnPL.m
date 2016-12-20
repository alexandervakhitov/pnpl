function [R, t, fail_flag, s] = OPnPL(U,u,xs, xe, Xs, Xe,label_polish)
    tm0 = tic;
%will polish the solution in default
    if nargin < 7
        label_polish = 'polish';
    end    
    n = size(U,2);
    
    %normalize lines
    Xd = Xe - Xs;
    if (size(Xd, 2) > 0)
        Xdn = sqrt(Xd(1, :).^2 + Xd(2, :).^2 + Xd(3, :).^2);
        Xd = Xd ./ repmat(Xdn, 3, 1);
        Xe = Xs + Xd;
    end

    l2d = opnp_lines.pnl_preprocess(xs, xe);

    tm1 = toc(tm0);
    tm0 = tic;
    [R0, t0, fail_flag1] = opnp_lines.opnpl_main1(U, u, Xs, Xe, l2d, label_polish);            
    tm2 = toc(tm0);
    tm0 = tic;
 
    [minInd Xsc Xec Lc] = opnp_lines.findBestRTReproj(R0, t0, xs, xe, Xs, Xe, U, u);
    [Xs1, Xe1] = opnp_lines.moveLineCloser2DLcVector(Xsc, Xec, xs, xe, Lc);
    R = R0(:, :, minInd);
    t = t0(:, minInd);    
    Xs = R'*Xs1 - repmat(R'*t, 1, size(Xs1, 2));
    Xe = R'*Xe1 - repmat(R'*t, 1, size(Xs1, 2));  

    tm3 = toc(tm0);
    tm0 = tic;
    [R0n, t0n, fail_flag2] = opnp_lines.opnpl_main1(U, u, Xs, Xe, l2d, label_polish);            
    
    nvar = size(R0, 3) + size(R0n, 3);
    R = zeros(3, 3, nvar);
    R(:, :, 1:size(R0, 3)) = R0;
    R(:, :, size(R0, 3)+1:end) = R0n;
    t = zeros(3, nvar);
    t(:, 1:size(t0, 2)) = t0;
    t(:, size(t0, 2)+1:end) = t0n;

    if (fail_flag1 == 1 || fail_flag2 == 1)
        fail_flag = 1;
    else
        fail_flag = 0;
    end
    tm4 = toc(tm0);
    s = [tm1, tm2, tm3, tm4];
end
function optInd = findBestRT(R01, R, t01, t)
    minDiff = 1e10;
    optInd = 1;
    for jj = 1:size(R01, 3)
        diff = norm(R01(:,:, jj) - R) + norm(t01(:, jj) - t);
        if (diff < minDiff)
            minDiff = diff;
            optInd = jj;
        end
    end
end

function [Xs_new, Xe_new] = moveLineCloser2D(Xsi, Xei, xs, xe)
    lineCoefs = cross(Xsi, Xei);
    lineCoefs = lineCoefs / norm(lineCoefs(1:2));    
    lineSpt = lineCoefs(1:2)*(-lineCoefs(3));
    lineDir = [-lineCoefs(2); lineCoefs(1)];
    ls = lineDir'*(xe-xs);
    if (ls<0)
        lineDir = -lineDir;
    end
    segLen = norm(xs-xe);
%     b = [xs-lineSpt; xe-lineSpt-lineDir*segLen];
%     A = [lineDir; lineDir];
%     sol = A\b;
    sol = lineDir'*(xs+xe-2*lineSpt)/2-segLen/2;
    xsn = lineSpt + sol*lineDir;
    xen = lineSpt + (sol+segLen)*lineDir;
    
    Xs_new = reproject2DtoLine(xsn, Xsi, Xei);
    Xe_new = reproject2DtoLine(xen, Xsi, Xei);
end

function Xs_new = reproject2DtoLine(xsi, Xsi, Xei)
    if (length(xsi) ~= size(xsi, 1))
        fprintf('ERROR in reproject2DtoLine \n');
    end
    if (length(xsi) == 2)
        xsi = [xsi; 1];
    end
    Xdir = Xei-Xsi;
    A = [xsi -Xdir];
    b = Xsi;
    sol = A\b;
    Xs_new = Xsi + Xdir*sol(2);
end

    
      
