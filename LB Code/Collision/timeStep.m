function [CL,LB] = timeStep(A,omega,dt,T,f,disc)

%Replace the value of tau
A = subs(A,'tau',T);
A = subs(A,'dt',dt);
omega = subs(omega,'tau',T);
omega = subs(omega,'dt',dt);
LB = matlabFunction(omega.*dt+f,'vars',[f]);
if (strcmp(disc,'explicit'))
%Step in time using an explicit Euler scheme
CL = A.*dt + eye(size(A));
else if (strcmp(disc,'implicit'))
CL = inv(-A.*dt + eye(size(A)));
end

end