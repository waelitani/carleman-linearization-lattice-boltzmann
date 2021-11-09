close all;
figure('DefaultAxesFontSize',14,'Name','Carlemann Linearization of Logistic Equation','WindowState','Maximized','DefaultLineLineWidth', 2);

colors = ['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30'];

subplot(2,Q,2*Q);
stopAt = 0; %fromend
stopOrder = length(order);
%stopOrder = 5;

%Set plotting limits
lwl = [0 0 0];
uwl = [1 1 1];

for q = 1:1:Q
subplot(2,Q,q);
hold on;
xlabel('Time (sec)');
ylabel(['f',num2str(q)]);



for o = 1:stopOrder
    if (order(o) < 10)
        st = ':';
        lw = 1.5;
        mm = 'x';
    else
        st = '-';
        lw = 1.5;
        mm = 'x';
    end
    plot((0:time-stopAt).*dt,VResults{o,whichT}(q,1:end-stopAt),'DisplayName',['Oc = ',num2str(order(o))],'LineStyle',st,'Marker',mm,'MarkerIndices',1:(1/dt):time,'LineWidth',lw,'Color',colors((mod(o,size(colors,1))>0)*(mod(o,size(colors,1)))+(mod(o,size(colors,1))==0)*(5),:));
    axis([0 time*dt lwl(q) uwl(q)]);
end
plot((0:time-stopAt).*dt,fResults{whichT}(q,1:end-stopAt),'DisplayName','Exact','Color','k');
%axis([0 time*dt 0 1]);
%yline(1, 'r--','DisplayName','Asymptote');
title(['f',num2str(q)]);
xticks(0:100:time);

xlabel('Time (sec)');

end

legend('Orientation','horizontal','Position',[0.1991 0.025 0.6 0.01]);

for q = 1:1:Q

	subplot(2,Q,Q+q);

    for o = 1:1:stopOrder
    if (order(o) < 10)
        st = ':';
        lw = 1;
        mm = 'none';
    else
        st = '-';
        lw = 1.5;
        mm = 'x';
    end
        er = abs(VResults{o,whichT}(q,1:end-stopAt)-fResults{whichT}(q,1:end-stopAt))./abs(fResults{whichT}(q,1:end-stopAt));
        plot((0:time-stopAt).*dt,er,'DisplayName',['Oc = ',num2str(order(o))],'LineStyle',st,'Marker',mm,'MarkerIndices',1:(1/dt):time,'LineWidth',lw,'Color',colors((mod(o,size(colors,1))>0)*(mod(o,size(colors,1)))+(mod(o,size(colors,1))==0)*(5),:));
            hold on;
    end

    %axis([0 time*dt 0 1]);
    if (q == 1)
        ylabel('Normalized Error');
    end
    axis([0 time*dt 0 0.1]);

end

sgtitle({'Carlemann Linearization of D1Q3',['Initial Velocity -> ',num2str(uu),' | dt/\tau ',num2str(eval(subs(1/T(whichT),dt,1)))]});

print([num2str(eval(subs(1/T(whichT),dt,1))),'.0.ps.png'],'-dpng');
