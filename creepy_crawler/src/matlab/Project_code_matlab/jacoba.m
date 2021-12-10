function J_a = jacoba(S,M,q)    
    % your code here 
    n = length(q);
    T_s = fkine(S,M,q, 'space');
    J_s = jacob0(S,q);
    J_a = J_s(4:6, 1:n) - skew(T_s(1:3,4)) * J_s(1:3, 1:n);    
end