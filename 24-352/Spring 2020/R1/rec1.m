%

clear all ; close all ; clc ; 

t = 0:0.1:5 ; 
x = (2/3)*exp(-5*t) + (1/3)*exp(4*t) ; 

figure(1)
plot(t,x,'-b') 
hold on

y0 = [1, 0] ; 
tspan = [0 5] ; 

m = 1 ; 
delB = -6 ; 
ch = 5 ; 

[t,y] = ode45(@(tspan,y0) odefun(tspan,y0,m,delB,ch),tspan,y0) ; 

figure(2)
plot(t,y(:,1)) 

figure(3) 
plot(t,delB*y(:,1),'-b')
hold on 
plot(t,ch*y(:,2),'-r')


%% 
clear all ; close all ; clc ; 
filename = '/home/emma/Downloads/01 - Love On Top.mp3' ; 
[y, fs] = audioread(filename) ; 

y = y(7837000:end,1) ; % this is where she starts to go up in pitch (about 90 sec from end)
dt = 1/fs ; 
t = 0:dt:(length(y)*dt)-dt ; 

figure(1)
plot(t,y) ; 
xlabel('sec') ; ylabel('amplitude') ; 

figure(2)
plot(psd(spectrum.periodogram,y,'Fs',fs,'NFFT',length(y))); % only plot to 1000Hz (human vocal range)
% axis([0 1000 -150 0]) ; 

%% Function definitions

function ydot = odefun(t,y,m,delB,ch)
    ydot = [y(2) ; (1/m) * ( (delB * y(1)) - (ch * y(2)) )] ; 
end