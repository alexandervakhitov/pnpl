function ret = inBox(Pts, Box)
    ret = 1;
    for i = 1:size(Pts, 2)
        for j = 1:3
            if (Pts(j, i) < Box(j, 1) || Pts(j, i) > Box(j, 2))
                ret = 0;
            end
        end
    end
end