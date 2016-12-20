function [err errl errt] = compute_reprojection(R1, t1, index_best, npt, nlines, XXw, xxn, Xs, Xe, xs, xe)
    reproj = R1(:,:,index_best)*XXw+t1(:,index_best)*ones(1,npt);
    reproj = reproj./repmat(reproj(3,:),3,1);
    errp = xxn-reproj(1:2,:);
    err = sqrt(sum(sum(errp.*errp))/npt);
    reproj_l =  reprojErrLines(R1(:,:,index_best), t1(:,index_best), Xs, Xe, xs, xe);
    errl = sqrt(sum(sum(reproj_l.*reproj_l))/2/nlines);
    errt = [errp reproj_l];
    errt = sqrt(sum(sum(errt.*errt))/(npt+nlines));
end