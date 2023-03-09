% Performance index to be minimized 
function J = PI(K,A,B,C,Q,R)
    Ac = A-B*K*C; 
    rmax =max(real(eig(Ac)));
    if rmax < 0
        P = lyap(Ac',Q+C'*K'*R*K*C);
        J = 0.5*trace(P);
    else
        J = 1e20 ;
    end
end