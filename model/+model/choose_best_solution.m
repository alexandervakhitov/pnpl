function [index_best, y] = choose_best_solution(R1, t1, R, t)
    index_best = 1;
    error = inf;
    y = 1e10*ones(2,1);
    for jjj = 1:size(R1,3)
        tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
        if sum(tempy) < error
            y = tempy;
            error = sum(tempy);
            index_best = jjj;
        end
    end
end