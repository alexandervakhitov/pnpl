function [R, t, XLTw, XXw, xxn, xs,xe,Xs,Xe] = setup_planar_scene(npt, nlines, segLenBig, segLenShiftScale, nll, nlp, fpix)


    XXw= [xrand(2,npt,[-2 2]); zeros(1,npt)];
    R= rodrigues(randn(3,1));
    t= [rand-0.5;rand-0.5;rand*4+4];
    Xc= R*XXw+repmat(t,1,npt);

    Xd = [xrand(1,nlines,[-2 2]); xrand(1,nlines,[-2 2]); zeros(1, nlines)];        


    XL1 = zeros(3, nlines);
    XL2 = zeros(3, nlines);
    for lineInd = 1:nlines
%             seglen = segLenPS(i)*randn;
        gend = 0;
        while (~gend)
            Xm = [xrand(1,1,[-2 2]); xrand(1,1,[-2 2]); 0];          
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
            segLen = (rand+0.5)*segLenBig;
            Xp1 = Xm2 + 0.5*Xdir*segLen;
            Xp2 = Xm2 - 0.5*Xdir*segLen;     
            XLT1(:, lineInd) = Xp1;
            XLT2(:, lineInd) = Xp2;

            Xp1 = R*Xp1 + t;
            Xp2 = R*Xp2 + t;

            XLp(1:2, lineInd) = Xp1(1:2)/Xp1(3)*fpix + randn(2,1)*nll;
            XLp(3:4, lineInd) = Xp2(1:2)/Xp2(3)*fpix + randn(2,1)*nll;
            if (inBox([Xp1 Xp2 R*XL1(:, lineInd)+t R*XL2(:, lineInd)+t], [-2 2; -2 2; 4 8]))
                gend = 1;
            end
        end
    end    

    XL = {XL1, XL2};

    XLw = {XL1, XL2};
    XLTw = {XLT1, XLT2};
    % projection        

    xx= fpix*[Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)];
    xxn= xx+randn(2,npt)*nlp;
    xxn = xxn / fpix;

    XLp = XLp / fpix; 

    xs = XLp(1:2, :);
    xe = XLp(3:4, :);
    Xs = XLw{1};
    Xe = XLw{2};
end