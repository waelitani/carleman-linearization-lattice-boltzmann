function [CL,LB] = timeStep(A,omega,dt,T,f)

%Replace the value of tau
A = subs(A,'tau',T);
omega = subs(omega,'tau',T);

%Step in time using an explicit Euler scheme
CL = A.*dt + eye(size(A,1));
LB = matlabFunction(omega.*dt+f,'vars',f);

end