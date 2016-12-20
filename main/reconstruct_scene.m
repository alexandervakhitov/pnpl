function reconstruct_scene(opts)
    for i = 1 : length(opts.fNamesBase)
        ims(:,:,:,i) = imread(opts.fNamesBase{i});
    end
    [gxs gxe] = getLineData3Frames(ims, opts.baseNames, opts.folder_to_save, opts.scale);
    
    for i = 1:length(opts.manualPointsFiles)
        fIn = fopen(opts.manualPointsFiles{i}, 'r');
        ptsRaw = fscanf(fIn, '%f '); 
        fclose(fIn);
        xx(:, :, i) = reshape(ptsRaw, 2, length(ptsRaw)/2);        
    end
    
    thr = opts.fmat_thr;
    [F, inliers] = ransacfitfundmatrix7(x1, x2, thr);
    projsMan = [xx(:, :, 1); xx(:, :, 2)];
    %visualize inliers
    visualizeMatchingN(zeros(4,0), zeros(4,0), projsMan, ims);

    E = K'*F*K;
    [R,b,votemin] = extractMotion(xx(:,:,1),xx(:,:,2),E,K,K);    
    rods = zeros(3,2);
    rods(:,2) = rodrigues(R);
    ts = zeros(3, 2);
    ts(:, 2) = b;    

    [rods, ts, pts3d, recPtInds] = runSbaPair(projsMan, rods, ts, opts.K, opts.distortion);
    
    filePtIn = fopen(opts.siftPointsFile, 'r');
    ptMatchRaw = fscanf(filePtIn, '%f ');
    fclose(filePtIn);
    ptMatchCrds = reshape(ptMatchRaw, 4, length(ptMatchRaw)/4)/scale;
    
    %triangulate pts
    R1 = eye(3);
    t1 = zeros(3,1);
    R2 = rodrigues(rods(:, 2));
    t2 = ts(:, 2);
    kc = zeros(10,1);
    kc(1) = opts.distortion;
    sift_p3ds = zeros(3, size(ptMatchCrds, 2));
    sift_gi = [];
    for pi = 1:size(ptMatchCrds, 2)
        proj1 = ptMatchCrds(1:2, pi);
        proj2 = ptMatchCrds(3:4, pi);
        [xn1,dxdf,dxdc,dxdk,dxdalpha] = normalize2(proj1,K(1,1)*ones(2,1),K(1:2, 3),kc,0);
        [xn2,dxdf,dxdc,dxdk,dxdalpha] = normalize2(proj2,K(1,1)*ones(2,1),K(1:2, 3),kc,0);
        [p3d, ret, err] = findPointPositionsLMNL([R1 t1], [R2 t2], xn1', xn2');
        sift_p3ds(:, pi) = p3d;
        sift_errs(pi) = err;
        if (err < opts.sift_reproj_err_thr)
            sift_gi = [sift_gi pi];
        end
    end
    
    %load and triangulate manual marked line segments
    for i = 1:length(opts.demoDetailsFiles)
        fIn = fopen(opts.demoDetailsFiles{i}, 'r');
        ptsRaw = fscanf(fIn, '%f '); 
        fclose(fIn);
        demo_xx(:, :, i) = reshape(ptsRaw, 2, length(ptsRaw)/2);        
    end
    R1 = rodrigues(rods(:,1));
    R2 = rodrigues(rods(:,2));
    t1 = ts(:, 1);
    t2 = ts(:, 2);
    for i = 1:size(demo_xx, 2)
        [p3d] = findPointPositionsLMNL(K*[R1 t1], K*[R2 t2], demo_xx(1:2, i, 1)', demo_xx(1:2, i, 2)');
        pts3d_demo(:, i) = p3d;
    end    
    
    
end