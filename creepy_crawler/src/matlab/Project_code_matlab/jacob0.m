function J = jacob0(S,q)
% Your code here
    J = [];
    T = eye(4);
    for i = 1:length(q)
        J = [J adjoint(S(:,i),T)];
        T = T * twist2ht(S(:,i), q(i));
    end
end