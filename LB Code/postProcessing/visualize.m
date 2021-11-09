%% Order Data
for i = 1:LO
    for j = 1:TO
        for k = 1:Kmax
            orderDat(i,j,k) = order(o);
        end
    end
end
    
%% Calculate distance from equilibrium

dist = zeros(TO,Kmax);

parfor k = 1:Kmax
	for j = 1:TO

            dist(j,k) = sqrt(sum((reshape(fResults{j,k}(:,1),[Q,1])-[1/6;4/6;1/6]).^2));
            
    end
end

dist = dist./max(max(dist));

%% Calculate engineered initial point
engData = zeros(LO,TO,Q);
enC = zeros(Q,1);
engDataC = zeros(LO,TO,Q);
enger = zeros(LO,TO);
engDataU = zeros(LO,TO);
engDataUC = zeros(LO,TO);
enguer = zeros(LO,TO);
engDataP = zeros(LO,TO);
engDataPC = zeros(LO,TO);
engper = zeros(LO,TO);

for i = 1:LO
    for j = 1:TO
    engData(i,j,:) = (LB{i,j}(1/6-uu/2,4/6,1/6+uu/2));
    
    enC = eval(CL{i,j}*subs(V{i},f,[1/6-uu/2;4/6;1/6+uu/2]));
    engDataC(i,j,:) = enC(end-Q:end-1);
    enger(i,j) = abs((engDataC(i,j)-engData(i,j))./engData(i,j));
    
    engDataU(i,j) = uF(engData(i,j,:));
    engDataUC(i,j) = uF(engDataC(i,j,:));
    enguer(i,j) = abs((engDataUC(i,j)-engDataU(i,j))./engDataU(i,j));
    
    engDataP(i,j) = rhoF(engData(i,j,:));
    engDataPC(i,j) = rhoF(engDataC(i,j,:));
    engper(i,j) = abs((engDataPC(i,j)-engDataP(i,j))./engDataP(i,j));
    end
end
    

%% Scatterplot
close all;
for j = 1:TO
figure('WindowState','maximized');
subplot(1,2,1);
title('Velocity');
hold on;
for o = 1:LO
scatter((o).*ones([size(uer,3),1]),uer(o,j,:),[],dist(j,:));
end

boxplot(reshape(uer(:,j,:),[Kmax,LO]),order,'orientation','vertical','MarkerStyle','none');

axis([0 LO+1 0 1]);
plot(enguer(:,j),'color','r','Linewidth',2,'DisplayName','Chosen Initial Condition');
colorbar();

subplot(1,2,2);
title('Density');
hold on;
for o = 1:LO
scatter((o).*ones([size(per,3),1]),per(o,j,:),[],dist(j,:));
end
boxplot(reshape(per(:,j,:),[Kmax,LO]),order,'orientation','vertical','MarkerStyle','none');
axis([0 LO+1 0 1]);
plot(engper(:,j),'color','r','Linewidth',2,'DisplayName','Chosen Initial Condition');
colorbar();

suptitle(['dt/\tau = ',num2str(oms(j)),'Initial Velocity ',num2str(uu),' -> | Random Sample Size = ',num2str(Kmax)]);

end