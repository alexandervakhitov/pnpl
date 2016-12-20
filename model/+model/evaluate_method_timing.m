function [R1, t1, is_fail, s] = evaluate_method_timing(meth, Xs, Xe, xs, xe, XXw, xxn, XLTw)
    is_fail = 0;
    s = 0;
    R1 = [];
    t1 = [];
    
    if (meth.met_type == 3 && meth.min_pt_num > size(XXw, 2))
        is_fail = 0;
        s = 0;
        return;
    end
    
    if (meth.met_type == 1 && meth.min_ln_num > size(Xs, 2))
        is_fail = 0;
        s = 0;
        return;
    end

    if (meth.met_type == 2 && meth.min_ln_num > size(Xs, 2)+size(XXw, 2))
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
            [R1,t1,ff,s] = meth.f(XXw, xxn, xs,xe,Xs,Xe);  
        end
        if (meth.met_type == 3) %points only                        
            [R1,t1] = meth.f(XXw, xxn);  
        end
        if (meth.met_type == 4) %points+true endpoints
            XXwf = [XXw XLTw{1} XLTw{2}];
            xxnf = [xxn xs xe];
            [R1,t1] = meth.f(XXwf, xxnf);  
        end
%         s = toc(time1);
    catch
        fprintf(['   The solver - ',meth.name,' - encounters internal errors! \n']);
        is_fail = 1;
%         index_fail = [index_fail, j];
%         break;
    end
end