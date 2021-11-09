function [Acoeffs, truncatedVars, truncatedNumVars] = carlemannLinearize(order,f,omega)
%Carlemann Linearization

%Choose the order of truncation
%order = 2;

vars = transpose(f.^([0:1:order]));

numVars = (order)+1;
C = sym(zeros(numVars,1));

%Expand the current matrix
%The ODE describing the new variables is obtained from:
%partial(var)/partial(t) = Sum_i partial(var)/partial(f_i) df_i/dt
parfor i = 1:1:numVars
    C(i,:) = [sum(gradient(vars(i),f).*omega)];
end

%Find the index of "truncated" variables up to the desired order
idx = find(polynomialDegree(vars,f) <= order);

%Determine the "truncated" variables
truncatedVars = vars(idx);

%and their corresponding differential equations
truncatedC =  C(idx,:);

clear idx C;

%Determine the number of truncated variables
truncatedNumVars = length(truncatedVars);

%Drop the higher order terms from the expressions of the chosen variables
parfor i = 1:1:truncatedNumVars
    [ms,vs] = coeffs(truncatedC(i),f,'All');
    vs = reshape(vs,[numel(vs),1]);
    ms = reshape(ms,[numel(ms),1]);
    poly = polynomialDegree(vs,f);
    idx = find(poly > order);
    ms(idx) = 0.*ms(idx);
    truncatedC(i) = sum(expand(ms.*vs));
end

clear i ms vs poly idx;

%Rearrange the system of equations in order of variables
tDegree = polynomialDegree(truncatedVars,f);
[tDegree,tSort] = sort(tDegree,'descend');

truncatedVars = truncatedVars(tSort);
truncatedC = truncatedC(tSort);

clear tDegree tSort;

%Replace the variables in terms of the discrete density
%with a new set of variables to obtain the coefficient matrix
%of the Carlemann linearization

%The last variable is not replaced since it corresponds to unity
%needed for constant terms of the corresponding polynomials
%Replace the variables with new variable set

%Declare new set of variables phi
phi = sym('phi',[truncatedNumVars-1,1]);

%Perform the symbolic substitution
for i = 1:truncatedNumVars-1
    parfor j = 1:1:numVars
    
    truncatedC(j) = subs(truncatedC(j),truncatedVars(i),phi(i));
    
    end
    
end

%Initialize the coefficients matrix
Acoeffs = sym(zeros(truncatedNumVars,truncatedNumVars));

for i = 1:1:truncatedNumVars-1

    [Acoeff,Avar] = coeffs(truncatedC(i),phi);
    for j = 1:1:length(Avar)
        idx = find([phi;1]==Avar(j));
        Acoeffs(i,idx) = Acoeff(j);
    end
    
end

clear i truncatedC Acoeff Avar j idx phi;
end