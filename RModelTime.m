close all; clear all; 
points= [2,8,20,50,200];
time=[1.063,9.314,145.83,2067.846,93321.886];

figure;
plot(points,time,'o')
xlabel('points'); ylabel('time (s)')
grid on 