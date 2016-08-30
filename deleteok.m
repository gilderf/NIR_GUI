
%figure(1)
clf;

x=[0:0.1:20]
y=sin(x)

subplot(1,1,1) ; % main plot

plot(x,y,'b')
xlabel('x','fontsize',18)
ylabel('sin(x)','fontsize',18)

axes('position',[0.2 0.7 0.2 0.1]) ; % 1st inset

plot(x,abs(y),'r')
xlabel('x','fontsize',15)
ylabel('abs(sin(x))','fontsize',15)

axes('position',[0.13 0.1 0.775 0.1]) ; % 2nd inset

plot(x,abs(y),'r')
%xlabel('x','fontsize',15)
%ylabel('abs(sin(x))','fontsize',15)