function AdT = adjoint(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
    AdT = [R zeros(3); skew(p)*R R];
end
