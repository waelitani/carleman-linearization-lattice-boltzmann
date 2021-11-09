%% Calculate Matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear the desk and focus on the job
close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lattice Boltzmann Formulation

%Determine the number of discrete velocities
Q = 3;
D = 1;

%Declare your weights
w = sym([1/6; 4/6; 1/6]);

%Work in lattice units
c = 1;
dt = 1;

%Declare numeric parameters
syms tau dt;

%Assume incompressible flow with unity constant density
rho = 1;

%Set your basis vectors for velocity determination
e = zeros(Q, Q);
e = [-1; 0; 1];

[f,omega,u,rho] = latticeBoltzmann(Q,D,w,e,dt,tau,c,rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up matrices
%Set tau for chosen omegas
T = 0:1:5;
T = 2-(10.^-T);
T = 1./T;
T = dt.*[Inf 10 5 2 T];
T = sort(T,'ascend');

%Set the number of timesteps to evolve the solution
time = 1;

%Set the number of iterations with randomized initial data
Kmax = 100;

%Choose your truncation order
%order = [1:5 10:5:25];
order = [1:5];
%Set the initial velocity
uu = 0.1;

%Choose initial conditions
%f0 = eval(w)+[-uu/2;0;+uu/2];

parfor o = 1:1:length(order)
    %display(['Order ',order(o)]);
    %Carlemann Linearization
    [A{o},V{o}] = carlemannLinearize(order(o),f,omega);
    
end
%% Replace Symbolics
dt = 1;
sizeOM = length(T);
parfor o = 1:1:length(order)
	for om = 1:sizeOM
        [CL{o,om},LB{o,om}] = timeStep(A{o},omega,dt,T(om),f,'explicit');
    end
end
%% Calculate Results
clc;
f0=0;
display('Computing Results');
[fResults,VResults] = computeResultsArray(CL,LB,V,order,f,f0,uu,omega,T,dt,time,Kmax);

%Post-process
dataTable;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%