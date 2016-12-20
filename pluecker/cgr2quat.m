function q = cgr2quat(s)


q = zeros(4,size(s,2));
for k = 1:size(s,2)
    q(:,k) = [s(:,k); 1];
    q(:,k) = q(:,k)/norm(q(:,k));
    if q(4,k) < 0
        q(:,k) = -q(:,k);
    end
end