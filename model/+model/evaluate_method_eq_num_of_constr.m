function [R1, t1, is_fail, s] = evaluate_method_eq_num_of_constr(meth, Xs, Xe, xs, xe, XXw, xxn, XLTw)
    is_fail = 0;
    s = 0;
    R1 = [];
    t1 = [];
    
    
     if (meth.min_pt_num > size(XXw, 2))
        is_fail = 0;
        s = 0;
        return;
    end
    
    if (meth.min_ln_num > size(Xs, 2))
        is_fail = 0;
        s = 0;
        return;
    end
    
    try
        time1 = tic;
        if (meth.met_type == 1) %lines only
            if (size(xs, 2) >= 3)
                [R1,t1] = meth.f(xs,xe,Xs,Xe);  
            else
                is_fail = 1;
                return;
            end
        end
        if (meth.met_type == 2) %lines+points only
            [XXw, xxn, Xs, Xe, xs, xe] = random_selection(XXw, xxn, Xs, Xe, xs, xe);
            [R1,t1] = meth.f(XXw, xxn, xs,xe,Xs,Xe);  
        end
        if (meth.met_type == 3) %points only                        
            [R1,t1] = meth.f(XXw, xxn);  
        end
        if (meth.met_type == 4) %points+true endpoints
           [XXw, xxn, XLTw{1}, XLTw{2}, xs, xe] = random_selection(XXw, xxn, XLTw{1}, XLTw{2}, xs, xe);
            XXwf = [XXw XLTw{1} XLTw{2}];
            xxnf = [xxn xs xe];
            [R1,t1] = meth.f(XXwf, xxnf);  
        end
        s = toc(time1);
    catch
        fprintf(['   The solver - ',meth.name,' - encounters internal errors! \n']);
        is_fail = 1;
%         index_fail = [index_fail, j];
%         break;
    end
end

function [XXw, xxn, Xs, Xe, xs, xe] = random_selection(XXw, xxn, Xs, Xe, xs, xe)
    pt_num = floor(size(XXw, 2) / 2);
    pts_inds = randsample(size(XXw, 2), pt_num);            
    XXw = XXw(:, pts_inds);
    xxn = xxn(:, pts_inds);
    lines_inds = randsample(size(XXw, 2), pt_num);
    Xs = Xs(:, lines_inds);
    Xe = Xe(:, lines_inds);
    xs = xs(:, lines_inds);
    xe = xe(:, lines_inds);
end