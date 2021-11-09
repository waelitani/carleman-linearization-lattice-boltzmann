%A = [-(3*dt)/(2*tau)	0	0	0; 0		-dt/tau	0	0; 0		0	-dt/(2*tau)	0; 0 0 0 0];
%B = [dt/(2*tau)	0	-dt/(4*tau)	dt/(6*tau);-dt/tau		0	-dt/tau		(2*dt)/(3*tau);-dt/(4*tau)	dt/(2*tau)	0	dt/(6*tau);0	0	0	0];
%F2 = abs(eval(eig(B)));
x = [];
for i = 0.5:0.1:10
    x = [x eval(subs(F2*tau,tau,i))];
end

figure;
hold on;

for i = 1:size(x,1)
    plot(0.5:0.1:10,2*x(i,:))
end

xlabel('\tau');
ylabel('||F2||/|Real(\lambda_1(F1))|');
title('Measure of Non-Linearity');