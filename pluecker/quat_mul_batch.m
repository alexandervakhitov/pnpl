function q = quat_mul_batch(q1,q2)



if (size(q1,2) == 1)
    q = zeros(4,size(q2,2));
    for k = 1:size(q2,2)
        q(:,k) = quat_mul(q1,q2(:,k));
    end
elseif (size(q2,2) == 1)
    q = zeros(4,size(q1,2));
    for k = 1:size(q1,2)
        q(:,k) = quat_mul(q1(:,k),q2);
    end
else
    assert(0, 'At least one of the quaternions should be one dimensional');
end
