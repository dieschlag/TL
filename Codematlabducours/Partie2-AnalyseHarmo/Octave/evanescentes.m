%%% S. Rossignol -- 05/08/13

clear;
close all;


%%% I

sigi1=[zeros(1,20) exp(-[-49.5:49.5].^2/400) zeros(1,180)];
timi=[0:length(sigi1)-1]-length(sigi1)/2+0.5;


%%% 1

f0=0.015;
sigi2=sigi1.*cos(2*pi*f0*timi);

figure(1);
clf;
plot(timi,sigi1,'Linewidth',3);
hold on;
plot(timi,sigi2,'Linewidth',3,'r');
title('onde evanescente (t=-80,f=0.015)','Fontsize',25);
ylim([-1.1 1.1]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;
print -depsc2 evanescentes1.eps


%%% 2

f0=0.030;
sigi2=sigi1.*cos(2*pi*f0*timi);

figure(2);
clf;
plot(timi,sigi1,'Linewidth',3);
hold on;
plot(timi,sigi2,'Linewidth',3,'r');
title('onde evanescente (t=-80,f=0.030)','Fontsize',25);
ylim([-1.1 1.1]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;
print -depsc2 evanescentes2.eps


%%% 3

f0=0.045;
sigi2=sigi1.*cos(2*pi*f0*timi);

figure(3);
clf;
plot(timi,sigi1,'Linewidth',3);
hold on;
plot(timi,sigi2,'Linewidth',3,'r');
title('onde evanescente (t=-80,f=0.045)','Fontsize',25);
ylim([-1.1 1.1]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;
print -depsc2 evanescentes3.eps


%%% 4

f0=0.060;
sigi2=sigi1.*cos(2*pi*f0*timi);

figure(4);
clf;
plot(timi,sigi1,'Linewidth',3);
hold on;
plot(timi,sigi2,'Linewidth',3,'r');
title('onde evanescente (t=-80,f=0.060)','Fontsize',25);
ylim([-1.1 1.1]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;
print -depsc2 evanescentes4.eps


%%% II

sigi1=[zeros(1,170) exp(-[-49.5:49.5].^2/400) zeros(1,30)];
timi=[0:length(sigi1)-1]-length(sigi1)/2+0.5;


%%% 1

f0=0.015;
sigi2=sigi1.*cos(2*pi*f0*timi);

figure(5);
clf;
plot(timi,sigi1,'Linewidth',3);
hold on;
plot(timi,sigi2,'Linewidth',3,'r');
title('onde evanescente (t=70,f=0.015)','Fontsize',25);
ylim([-1.1 1.1]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;
print -depsc2 evanescentes5.eps


%%% 2

f0=0.030;
sigi2=sigi1.*cos(2*pi*f0*timi);

figure(6);
clf;
plot(timi,sigi1,'Linewidth',3);
hold on;
plot(timi,sigi2,'Linewidth',3,'r');
title('onde evanescente (t=70,f=0.030)','Fontsize',25);
ylim([-1.1 1.1]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;
print -depsc2 evanescentes6.eps


%%% 3

f0=0.045;
sigi2=sigi1.*cos(2*pi*f0*timi);

figure(7);
clf;
plot(timi,sigi1,'Linewidth',3);
hold on;
plot(timi,sigi2,'Linewidth',3,'r');
title('onde evanescente (t=70,f=0.045)','Fontsize',25);
ylim([-1.1 1.1]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;
print -depsc2 evanescentes7.eps


%%% 4

f0=0.060;
sigi2=sigi1.*cos(2*pi*f0*timi);

figure(8);
clf;
plot(timi,sigi1,'Linewidth',3);
hold on;
plot(timi,sigi2,'Linewidth',3,'r');
title('onde evanescente (t=70,f=0.060)','Fontsize',25);
ylim([-1.1 1.1]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;
print -depsc2 evanescentes8.eps

