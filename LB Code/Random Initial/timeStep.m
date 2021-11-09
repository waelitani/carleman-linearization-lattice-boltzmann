function [CL,LB] = timeStep(A,omega,dt,T,f,disc)

%Replace the value of tau
A = subs(A,'tau',T);
omega = subs(omega,'tau',T);

if (strcmp(disc,'explicit'))
%Step in time using an explicit Euler scheme
CL = A.*dt + eye(size(A));
LB = matlabFunction(omega.*dt+f,'vars',f);
else if (strcmp(disc,'implicit'))
CL = A.*dt + eye(size(A));
LB = matlabFunction(omega.*dt+f,'vars',f);
end

end