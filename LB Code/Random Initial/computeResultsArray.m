function [fResults,VResults] = computeResultsArray(CL,LB,V,order,f,f0,uu,omega,T,dt,time,Kmax)

%Recompute number of discrete velocities
Q = length(f);

for o = 1:1:length(order)
    
    for om = 1:1:length(T)

    parfor k = 1:1:Kmax
		
%  			fL = rand*(1-uu);
%             fR = uu+fL;
%             fC = 1-fR-fL;
%  			f0 = [fL; fC; fR];
            f0 = rand(3,1);
            f0 = f0./sum(f0);

            %Store the initial value in the results array for the respective tau
            fT = sym(zeros(length(f),time+1));
			fT(:,1) = f0;
			V0 = subs(V{o},f,f0);
			VT = sym(zeros(length(V{o}),time+1));
			VT(:,1) = V0;

            if (o == 1)
                %Step in time
                for j = 2:1:time+1

                  fT(:,j) = LB{o,om}(fT(1,j-1),fT(2,j-1),fT(3,j-1));
                  
                end
                
                %Store the results in a cell array
                fResults{om,k} = fT;

            end
            
            for j = 2:1:time+1

                %VT(:,j) = (CL{o,om}^(j-1))*V0;
				VT(:,j) = CL{o,om}*VT(:,j-1);
            end
            
            VResults{o,om,k} = VT(end-Q:end-1,:);

        end

    end
end

end