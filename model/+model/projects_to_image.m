function is_ok = projects_to_image(projs, resx, resy)
    is_ok = 1;
    for i = 1:size(projs, 2)
        if (abs(projs(1, i)) > resx/2 || abs(projs(2, i)) > resy/2)
            is_ok = 0;
        end
    end
end