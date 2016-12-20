function [R, t, XLTw, XXw, xxn, xs,xe,Xs,Xe] = setup_3d_scene(npt, nlines, segLenBig, segLenShiftScale, nll, nlp, fpix, is_len_var)        
        if (nargin < 8)
            is_len_var = 1;
        end
        
        res_x = 640;
        res_y = 480;
        % generate 3d coordinates in camera space
        
        R= rodrigues(randn(3,1));
        
        pts_ok = 0;
        while (~pts_ok)
            Xc = [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];                              
            xx= fpix*[Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)];                           
            xxn= xx+randn(2,npt)*nlp;                        
            pts_ok = model.projects_to_image(xxn, res_x, res_y);
        end
        
        xxn = xxn / fpix;            
        
        Xd = [xrand(1,nlines,[-2 2]); xrand(1,nlines,[-2 2]); xrand(1,nlines,[-2 2])];        
        
        XL1 = zeros(3, nlines);
        XL2 = zeros(3, nlines);
        XLT1 = zeros(3, nlines);
        XLT2 = zeros(3, nlines);
        XLp = zeros(4, nlines);
        

        
        for lineInd = 1:nlines
%             seglen = segLenPS(i)*randn;
            gend = 0;
            while (~gend)
                Xm = [xrand(1,1,[-2 2]); xrand(1,1,[-2 2]); xrand(1,1,[4 8])];          
                Xdir = Xd(:, lineInd);
                Xdir = Xdir / norm(Xdir);                                

                segLen = (rand+0.5)*segLenBig;%*0.75*(rand+1); 
                segShift = segLenShiftScale*(2*rand-1)*0.5*Xdir;
                Xm1 = Xm + segShift;
                Xp1 = Xm1 + 0.5*Xdir*segLen;
                Xp2 = Xm1 - 0.5*Xdir*segLen;              
                XL1(:, lineInd) = Xp1;
                XL2(:, lineInd) = Xp2;                    
    
                segShift = segLenShiftScale*(2*rand-1)*0.5*Xdir;
                Xm2 = Xm + segShift;
                if (is_len_var)
                    segLen = (rand+0.5)*segLenBig;
                end
                Xp1 = Xm2 + 0.5*Xdir*segLen;
                Xp2 = Xm2 - 0.5*Xdir*segLen;     
                XLT1(:, lineInd) = Xp1;
                XLT2(:, lineInd) = Xp2;

                XLp(1:2, lineInd) = Xp1(1:2)/Xp1(3)*fpix + randn(2,1)*nll;
                XLp(3:4, lineInd) = Xp2(1:2)/Xp2(3)*fpix + randn(2,1)*nll;                                
                if (inBox([Xp1 Xp2 XL1(:, lineInd) XL2(:, lineInd)], [-2 2; -2 2; 4 8]) && ...
                        model.projects_to_image([XLp(1:2, :) XLp(3:4, :)], res_x, res_y))
                    gend = 1;
                end                
            end
        end    

        XL = {XL1, XL2};
        t = mean([Xc XL1 XL2], 2);
        XXw= inv(R)*(Xc-repmat(t,1,npt));
        XLw = {inv(R)*(XL1-repmat(t,1,nlines)), inv(R)*(XL2-repmat(t,1,nlines))};
        XLTw = {inv(R)*(XLT1-repmat(t,1,nlines)), inv(R)*(XLT2-repmat(t,1,nlines))};


        
%         XLp = {[XL1(1,:)./XL1(3,:); XL1(2,:)./XL1(3,:)]*fpix + randn(2,nlines)*nl, [XL2(1,:)./XL2(3,:); XL2(2,:)./XL2(3,:)]*fpix + randn(2,nlines)*nl};


        
        XLp = XLp / fpix; 
        
        xs = XLp(1:2, :);
        xe = XLp(3:4, :);
        Xs = XLw{1};
        Xe = XLw{2};
end