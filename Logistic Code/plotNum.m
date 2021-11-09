figure('DefaultAxesFontSize',14,'Name','Carlemann Linearization of Logistic Equation','WindowState','Maximized','DefaultLineLineWidth', 2);

hold on;
xlabel('Time (sec)');
ylabel('x');

plot((0:time).*dt,fResults{1},'DisplayName','Exact','Color','k');

for o = 1:length(order)
    if (order(o) < 10)
        st = ':';
        lw = 1.25;
        mm = 'none';
    else
        st = '-';
        lw = 1.5;
        mm = 'x';
    end
    plot((0:time).*dt,VResults{o,1},'DisplayName',['N = ',num2str(order(o))],'LineStyle',st,'Marker',mm,'MarkerIndices',1:(1/dt):time,'LineWidth',lw);
end

axis([0 time*dt 0 2]);
yline(1, 'r--','DisplayName','Asymptote');
title('Numerical Solution');
xticks(0:1:5);

legend('Orientation','horizontal','Position',[0.2177 0.865 0.6 0.05]);
suptitle({'Carlemann Linearization of Logistic Equation',['Initial condition x(t = 0) = ',num2str(f0),' | Timestep ',num2str(dt), ' sec']});

print([num2str(f0),'.0.png'],'-dpng');
