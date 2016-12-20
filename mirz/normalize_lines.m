function lines = normalize_lines(lines,fc,cc,alpha_c)

% $Id: normalize_lines.m 1373 2012-06-10 06:29:39Z faraz $

N_ln = length(lines);

for k = 1:N_ln
    
    if size(lines(1).point1,1) == 1
        x1 = normalize(lines(k).point1',fc,cc,0,alpha_c); 
        x2 = normalize(lines(k).point2',fc,cc,0,alpha_c);
    else
        x1 = normalize(lines(k).point1,fc,cc,0,alpha_c); 
        x2 = normalize(lines(k).point2,fc,cc,0,alpha_c);
    end
    
    lines(k).pnt1_n = x1;
    lines(k).pnt2_n = x2; 
    
    M = [x1' 1;
        x2' 1];
    
    lines(k).nmoment = null(M); %normalized (measured) moment
    %lines(k).moment = nn/(norm(nn(1:2)));
    lines(k).no = k;
end
