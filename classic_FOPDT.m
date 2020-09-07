clc
clear all
close all

load('dados.mat')
to = length(y)/(max(t));
T=t(2)-t(1);
K=max(y)-min(y);

%% plot inicial

figure
plot(t, y)
hold on
plot(t,u)
hold off
grid
legend('Saída: y(t)', 'Entrada: u(t)');
xlabel('tempo(s)');
ylabel('Amplitude')
axis([1 20 3 62]);


%% Ziegler & Nichols

d1y = gradient(y,t);                                            % Numerical Derivative
d2y = gradient(d1y,t);                                          % Numerical Second Derivative
t_infl = interp1(d1y(21:end), t(21:end), max(d1y));             % Find ‘t’ At Maximum Of First Derivative
y_infl = interp1(t(21:end), y(21:end), t_infl);                 % Find ‘y’ At Maximum Of First Derivative
slope  = interp1(t, d1y, t_infl);                               % Slope Defined Here As Maximum Of First Derivative
intcpt = y_infl - slope*t_infl;                                 % Calculate Intercept
tngt = slope*t + intcpt;                                        % Calculate Tangent Line
figure(1)
plot(t, y,'LineWidth',1)
hold on
plot(t, tngt, '-r')                                             % Plot Tangent Line
hold on
plot(t_infl, y_infl, 'blackx')                                  % Plot Maximum Slope

p = polyfit(t, tngt, 1);                                        %acha a equacao da reta tangente

tau_max=(60-p(2))/p(1);                                         %valor extremo do tau
theta=(10-p(2))/p(1);                                           %valor do theta
tau=tau_max-theta;                                              %tau = tau_max-theta

% aqui criamos nossa transfer function em funcao dos valores encontrados e
% plotamos
s=tf('s');
num=K*exp(-theta*s);
den=tau*s+1;
sys=num/den;
sysD = c2d(sys,T,'zoh');

opt = stepDataOptions('InputOffset',0.2);
[y_ziegler,t_ziegler] = step(sysD, 20,opt);
plot(t_ziegler,y_ziegler, 'k.')
hold off

grid
legend('y(t)', 'Tangente','ponto inflexão', 'ziegler &  nichols')
axis([1 8 8 62])
mape = mean((abs(y_ziegler-y))./y)*100
title(['Ziegler & Nichols'])
xlabel('tempo(s)');
ylabel('Amplitude')

%% smith

y2 = 50*0.632+10; %ponto 1 em 64,2%
[ d, ix ] = min( abs( y-y2 ) );
y2 = y(ix);
t2 = t(ix);
y1 = 50*0.283+10; %ponto 1 em 28,3%
[ d, ix ] = min( abs( y-y1 ) );
y1 = y(ix);
t1 = t(ix);

tau_smith = (3/2)*(t2-t1);
theta_smith = t2-tau_smith;

figure(2)
plot(t, y,'LineWidth',1)
hold on
plot([t1 t2],[y1 y2], 'x')

% aqui criamos nossa transfer function em funcao dos valores encontrados e
% plotamos
s=tf('s');
num=K*exp(-theta_smith*s);
den=tau_smith*s+1;
sys=num/den;
sysD = c2d(sys,T,'zoh');

opt = stepDataOptions('InputOffset',0.2);
[y_smith,t_smith] = step(sysD, 20,opt);
plot(t_smith,y_smith, 'k.')
hold off

grid
mape = mean((abs(y_smith-y))./y)*100
legend('y(t)', 'coordenadas smith', 'Curva smith')
axis([1 8 8 62])
title('Método de smith')
xlabel('tempo(s)');
ylabel('Amplitude')

%% Sundaresan e Krishnaswamy

y2 = 50*0.853+10; %ponto 1 em 85,3%
[ d, ix ] = min( abs( y-y2 ) );
y2 = y(ix);
t2 = t(ix);
y1 = 50*0.353+10; %ponto 1 em 35,3%
[ d, ix ] = min( abs( y-y1 ) );
y1 = y(ix);
t1 = t(ix);

tau_sund = 0.67*(t2-t1);
theta_sund = 1.3*t1-0.29*t2;

figure(3)
plot(t, y,'LineWidth',1)
hold on
plot([t1 t2],[y1 y2], 'x')

% aqui criamos nossa transfer function em funcao dos valores encontrados e
% plotamos
s=tf('s');
num=K*exp(-theta_sund*s);
den=tau_sund*s+1;
sys=num/den;
sysD = c2d(sys,T,'zoh');

opt = stepDataOptions('InputOffset',0.2);
[y_sund,t_sund] = step(sysD, 20,opt);
plot(t_sund,y_sund, 'k.')
hold off

grid
mape = mean((abs(y_sund-y))./y)*100
legend('y(t)', 'coordenadas Sundaresan', 'Curva Sundaresan')
axis([1 8 8 62])
title('Método de Sundaresan e Krishnaswamy')
xlabel('tempo(s)');
ylabel('Amplitude')



%% Nishikawa

A0=0;
for i =1:length(t)
    A0=(K-(y(i)-y(1)))*T+A0;
end
T1=(A0/K)/T;

A1=0;
for i =1:T1
    A1=(y(i)-y(1))*T+A1;
end

tau_nish = (A1/(0.368*K));
theta_nish = (T1*T-tau_nish);

figure(4)
plot(t, y,'LineWidth',1)
hold on

% aqui criamos nossa transfer function em funcao dos valores encontrados e
% plotamos
s=tf('s');
num=K*exp(-theta_nish*s);
den=tau_nish*s+1;
sys=num/den;
sysD = c2d(sys,T,'zoh');

opt = stepDataOptions('InputOffset',0.2);
[y_nick,t_nick] = step(sysD, 20,opt);
plot(t_nick,y_nick, 'k.')
hold off

grid
mape = mean((abs(y_nick-y))./y)*100
legend('y(t)', 'Curva Nishikawa')
axis([1 8 8 62])
title('Método de Nishikawa')
xlabel('tempo(s)');
ylabel('Amplitude')


%% TODOS
figure(5);
plot(t,y,'LineWidth',1)
hold on
plot(t_ziegler,y_ziegler,'LineWidth',1)
plot(t_smith,y_smith,'LineWidth',1)
plot(t_sund,y_sund,'LineWidth',1)
plot(t_nick,y_nick,'LineWidth',1)
axis([1 8 8 62])
title('Métodos juntos')
xlabel('tempo(s)');
ylabel('Amplitude')
grid('on')
legend('original', 'Ziegler', 'Smith', 'Sundaresan', 'Nishikawa')