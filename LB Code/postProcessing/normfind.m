F1norm = 0;
F2norm = 0;
omega1 = subs(omega1,[dt tau],[1 1]);
omega2 = subs(omega2,[dt tau],[1 1]);

for q = 1:Q
    F1norm = max(F1norm,max(abs(coeffs(omega1(q),f))));
    F2norm = max(F2norm,max(abs(coeffs(omega2(q),f))));
end
    