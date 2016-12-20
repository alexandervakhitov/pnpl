function [R2,t2] = calcampose(XXc,XXw)
% XXc: 3D coordinates in camera frame. Its size is 3 by n.
% XXw: 3D coordinates in world frame. Its size is 3 by n.
% R2, t2: the pose of the object in the camera frame.

n= length(XXc);

X= XXw;%B
Y= XXc;%A

K= eye(n)-ones(n,n)/n;

ux= mean(X,2);
uy= mean(Y,2);
sigmx2= mean(sum((X*K).^2));

SXY= Y*K*(X')/n;
[U, D, V]= svd(SXY);
S= eye(3);
if det(SXY) < 0
    S(3,3)= -1;
end

R2= U*S*(V');
c2= trace(D*S)/sigmx2;
t2= uy-c2*R2*ux;

X= R2(:,1);
Y= R2(:,2);
Z= R2(:,3);
if norm(xcross(X,Y)-Z) > 2e-2
    R2(:,3)= -Z;
end

return

function c = xcross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
 
 return
