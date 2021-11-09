function [YSol,Y0] = solveA(A,V,f0)

    f = symvar(V);
    Y0 = eval(subs(V,f,f0));
    syms("Y(t)",[length(V) 1]);
    odes = diff(Y) == A*Y;
    YSol = dsolve(odes,Y(0) == Y0);
    YSol = YSol.(['Y',num2str(length(V)-1)]);

end