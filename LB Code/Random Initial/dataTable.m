%% Compute
%dooe = zeros(Kmax*length(T)*length(order)*time,4+Q+1);
dooe = zeros(1,4+Q+1);
oms = eval(1./subs(T,dt,1));

for k = 1:Kmax
    for om = 1:length(T)-1
        dist = sqrt(sum((fResults{om,k}(:,1)-eval(w)).^2));  
        for o = 2:length(order)
            for t = 2:time+1
                %time distance order omega error
                eqq = reshape(abs((VResults{o,om,k}(:,t)-fResults{om,k}(:,1))./fResults{om,k}(:,1)),[1 Q]);
                eqs = sqrt(sum(eqq.^2));
                dooe = [dooe;t-1 dist order(o) oms(om) eqq eqs];
            end
        end
    end
end
dooe(isinf(dooe)|isnan(dooe)) = 0;
%% Plot
close;
colormap(parula(16));
scatter3(dooe(:,2),dooe(:,3),dooe(:,4),40,log(dooe(:,4+Q+1)),'filled');
view(-31,14);
xlabel('Distance');
ylabel('Order');
yticks([2:5]);
zlabel('dt/\tau');
cb = colorbar;
cb.Label.String = 'RSS Normalized Error (log)';
title('Carlemann Linearization Error for D1Q3 Lattice Boltzmann');
caxis([-3 3]);
