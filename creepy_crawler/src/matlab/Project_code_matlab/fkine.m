function T = fkine(S,M,q,frame)
    % your code here
    if nargin < 4
        frame = 'space';
    end
    temp = eye(4);
    for i = 1:length(q)
        temp = temp * twist2ht(S(:,i),q(i));
    end
    if frame == "space"
        T = temp*M;
    end
    if frame == "body"
        T = M*temp;
    end
    % If needed, you can convert twists to homogeneous transformation matrices with:
    % twist2ht(S(i),q(i));
end