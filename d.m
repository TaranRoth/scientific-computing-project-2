clearvars
data = load("Hare_lynx");
t0 = 1907;
tf = 1935;
IC = [25.37 6.30];
params = [0.4;0.02;0.6;0.02];
costfixed = @(v) cost(v, t0, tf, data);
options = optimset('MaxFunEvals', 20000, 'MaxIter', 10000);
[final_params, err] = fminsearch(costfixed, [(IC');params], options);

[tvals, yvals] = ode45(@(t,y) f(t,y,final_params(3:6)),[t0,tf],final_params(1:2)');
plot(tvals,yvals(:,1),'Color','r')
hold on
plot(tvals,yvals(:,2),'Color','b')
legend('Hares (fit)','Lynxes (fit)', 'Hares (data)', 'Lynxes (data)')
xlabel("Time (years)")
ylabel("Population (thousands)")
hold on
xdata=data(62:90,1);
hdata=data(62:90,2);
ldata=data(62:90,3);
plot(xdata,hdata,'Marker','o','Color','r','MarkerFaceColor','r','LineStyle','none')
hold on
plot(xdata,ldata,'Marker','o','Color','b','MarkerFaceColor','b','LineStyle','none')
hold off

function result = f(t,y,p)
    h = y(1); l = y(2);
    a = p(1); b = p(2); c = p(3); d = p(4);
    result = [a*h - b*h*l; -c*l + d*h*l];
end

function y=cost(v, t0, tf, data)
    [tvals, yvals] = ode45(@(t,y) f(t,y,v(3:6)),[t0,tf],v(1:2)');
    sum = 0;
    for year=t0:tf
        year_idx = find(data == year);
        [min_val, idx] = min(abs(tvals - year));
        sum = sum + norm(yvals(idx, :) - data(year_idx,2:3));
    end
    y=sum;
end
