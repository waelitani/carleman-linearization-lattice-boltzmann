function [fResults,VResults] = computeResultsArray(A,V,order,f,f0,omega,T,dt,time,Kmax)

%Recompute number of discrete velocities
Q = length(f);
%timeO = time;
for o = 1:1:length(order)
    
    for om = 1:1:length(T)
        %time = timeO/min(1,T(om));
        %Construct MATLAB functions of lattice Boltzmann and the corresponding
        %Carlemann linearization for the chosen tau independent
        %of symbolic calculation to reduce computational cost
        [CL,LB] = timeStep(A{o},omega,dt,T(om),f);

        for k = 1:1:Kmax

            %Store the initial value in the results array for the respective tau
            fT = zeros(Q,time+1);
            fT(:,1) = f0;
            VT = (zeros(length(V{o}),time+1));
            V0 = subs(V{o},f,f0);
            VT(:,1) = V0;

            if (o == 1)
                %Step in time
                for j = 2:1:time+1

                  fT(:,j) = [LB(fT(j-1))];
                  
                end
                
                %Store the results in a cell array
                fResults{om,k} = fT;

            end
            
            parfor j = 2:1:time+1

                VT(:,j) = [(CL^(j-1))*V0];

            end
            
            VResults{o,om,k} = VT(end-Q:end-1,:);

        end

    end
end

end