close all;

Q = 1;

numRows = 2;
numCols = Q;

titleQ = {'Population'};

stopAt = 0;
stopOrder = length(order);

for om = 1:1:length(T)
    
    figure;
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);

    for q = 1:1:Q
    
        subplot(numRows,numCols,(numRows-1)*numCols+q);
        
        if (q == 1)
            ylabel('Normalized Error');
        end

        hold on;

        for o = 1:1:stopOrder
            
            semilogy(0:time-stopAt,abs(VResults{o,om}(q,1:end-stopAt)-fResults{om}(q,1:end-stopAt))./abs(fResults{om}(q,1:end-stopAt)));

        end
        
        axis([0 10 0 1]);

    end

    axis([0 10 0 1]);
    
    titleStr = ['Initial Population ', num2str(fResults{om}(:,1)), ' Omega ', num2str(dt*T(om)),' dt ',num2str(dt)];
    suptitle(titleStr);
    
    print([num2str(f0),' ',num2str(T(om)),'.0.pp.png'],'-dpng');
    close all;
    
end