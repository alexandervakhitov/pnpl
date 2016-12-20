function reproj_l = reprojErrLines(R, t, Xs, Xe, xs, xe)
    nlines = size(Xs, 2);
    reproj_l1 = R*Xs+t*ones(1,nlines);
    reproj_l1 = reproj_l1 ./ repmat(reproj_l1(3,:), 3, 1);
    reproj_l2 = R*Xe+t*ones(1,nlines);
    reproj_l2 = reproj_l2 ./ repmat(reproj_l2(3,:), 3, 1);

    lineEqs = [];
    for lineInd = 1:nlines
        c = cross(reproj_l1(:, lineInd), reproj_l2(:, lineInd));
        c = c / norm(c(1:2));
        lineEqs = [lineEqs c];
    end

    reproj_l = zeros(2, nlines);
    for lineInd = 1:nlines
        reproj_l(1, lineInd) = lineEqs(:, lineInd)'*[xs(:, lineInd); 1];
        reproj_l(2, lineInd) = lineEqs(:, lineInd)'*[xe(:, lineInd); 1];                    
    end
end