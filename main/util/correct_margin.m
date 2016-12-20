function correct_margin()
    pos = get(gca, 'Position');
    pos(2) = 0.2;
    pos(4) = 0.6;
    set(gca, 'Position', pos)
end