function xout = NextDeg2(xin)

% $Id: NextDeg2.m 1264 2012-01-20 02:03:57Z faraz $

len = length(xin);

x = xin(len:-1:1);

for k = len-1:-1:1

    if x(k) > 0
        
        x(k) = x(k) - 1;
        x(k+1) = sum(x(k+1:len)) + 1;
        x(k+2:len) = 0;
        
        xout = x(len:-1:1);
        
        return;
        
    end
end

x(1) = sum(x) + 1;
x(2:len) = 0;
xout = x(len:-1:1);


% len = length(xin);
% 
% x = xin(len:-1:1);
% if x(4) > 0
%     x(4) = x(4) - 1;
%     x(5) = x(5) + 1;
%     
% elseif x(3) > 0
%     x(3) = x(3) - 1;
%     x(4) = x(4) + x(5) + 1;
%     x(5) = 0;
%     
% elseif x(2) > 0
%     x(2) = x(2) - 1;
%     x(3) = x(3) + x(4) + x(5) + 1;
%     x(4) = 0;
%     x(5) = 0;
%     
% elseif x(1) > 0
%     x(1) = x(1) - 1;
%     x(2) = x(2) + x(3) + x(4) + x(5) + 1;
%     x(3) = 0;
%     x(4) = 0;
%     x(5) = 0;
%     
% else
%     x(1) = sum(x) + 1;
%     x(2:5) = 0;
%     
% end
% xout = x(len:-1:1);
