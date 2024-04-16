function out = vector_angles_R02(v1,v2)
% Takes in Nx3 matrices giving vector elements for v1 and v2. Each row of
% v1 and v2 is a new vector.
%out = acos((v1(:,1).*v2(:,1) + v1(:,2).*v2(:,2) + v1(:,3).*v2(:,3))./dot(vecnorm(v1')',vecnorm(v2')',2));

num = v1(:,1).*v2(:,1) + v1(:,2).*v2(:,2) + v1(:,3).*v2(:,3);
obtuse_check = num < 0;
den = (vecnorm(v1')'.*vecnorm(v2')');
out = acos(num./den);
out(obtuse_check) = out(obtuse_check) + pi/4;
end