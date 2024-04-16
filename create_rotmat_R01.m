function R = create_rotmat_R01(v1,v2)
    % Generates a rotation matrix to rotate from vector 1 to vector 2.
    v1 = v1/norm(v1); % Normalize
    v2 = v2/norm(v2); % Normalize
    e = cross(v1,v2); % Create rotation vector
    if sum(boolean(e))==0
        fprintf('Parallel vectors. Input intermediate rotation to clarify.\n')
        return
    end
    e = e/norm(e); % Normalize
    angle = acos(dot(v1,v2)); % Find angle
    c = cos(angle); s = sin(angle);
    % Rodriques' rotation formula
    K = [0,-e(3),e(2);e(3),0,-e(1);-e(2),e(1),0];
    R = eye(3)+s*K+(1-c)*(K^2);
end