%% Set up matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear the desk and focus on the job
close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lattice Boltzmann Formulation

dt = 0.1;

%Declare numeric parameters
syms tau;

syms f t;

omega = f*(1-f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set tau for chosen omegas
T = [0.1 0.2 0.5 1 1.1:0.1:5 10];

%Set the number of timesteps to evolve the solution
time = 5/dt;

%Set the number of iterations with randomized initial data
Kmax = 1;

%Choose your truncation order
order = [1:5 10:5:25 50];

%Carlemann Linearization    
for o = 1:1:length(order)
	
	order(o)
	[A{o},V{o}] = carlemannLinearize(order(o),f,omega);
	
end
%% Solve
LO = length(order);
parfor initCond = 0:10   

	%Choose initial conditions
	f0 = initCond/10;

	[fR{initCond+1},VR{initCond+1}] = computeResultsArray(A,V,order,f,f0,omega,T,dt,time,Kmax);
    
    for o = 1:1:LO
        [analyticSol{o,initCond+1}] = solveA(A{o},V{o},f0);
    end

end
%% Plot
    for initCond = 0:1:10   

	%Choose initial conditions
	f0 = initCond/10;
    if initCond ~= 0
    Tin = log(1+(1/f0));
    end
    fResults = fR{initCond+1};
    VResults = VR{initCond+1};
    plotSolReg;
    close all;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%