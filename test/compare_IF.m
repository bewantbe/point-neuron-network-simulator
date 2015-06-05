
x = load('neu_state_init.txt');
y = load('volt.txt');

figure(1);
m=1:50;
plot(m, x(m,1), m, y(m,1))

figure(2);
m=1:50;
plot(m, x(m,4), m, y(m,4))

figure(11);
m=1:1000;
plot(m, x(m,1)-y(m,1));
figure(12);
m=1:1000;
plot(m, x(m,4)-y(m,4));



figure(3);
m=1:50;
plot(m, x(m,2), m, y(m,2))

plot(m, x(m+1,2)-0.77880*x(m,2), m, y(m+1,2)-0.77880*y(m,2))


figure(4);
plot(m, x(m+1,5)-0.77880*x(m,5), m, y(m+1,5)-0.77880*y(m,5))

