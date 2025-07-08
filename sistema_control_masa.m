clc
close all
clear
 
%Sistema de Control de un sistema MASA FRICCION RESORTE
%datos
M1=5;% Kg
M2=4;% Kg
B1=10;  B2=10;  Ba=12;
K1=20;  K2=20;  Ka=16;
 
%% Representaciom mediante variables de estados
%x1=y1; x2=y2; x3=y'1; x4=y'2
%y=y1;  u=f1
A=[  0             0          1         0
     0             0          0         1
  -(K1+Ka)/M1   Ka/M1  -(B1+Ba)/M1    Ba/M1
     Ka/M2  -(Ka+K2)/M2    Ba/M2   -(Ba+B2)/M2]
 B=[0;0;1/M1;0]
 C=[0 1 0 0]
 D=0
 
 %% Sinulacion en lazo abierto tiempo continuo
 % Ante un escalon
 t=0:0.0001:6;
 [Yt,Xt]=step(A,B,C,D,1,t);
 figure
 plot(t,Yt)
 figure
 plot(t,Xt)
 legend('x1(t)','x2(t)','x3(t)','x4(t)')
 
 %% Discretizacion con Ts=0.1 seg y retenedor de orden cero
 Ts=0.1;
 [G,H,C,D]=c2dm(A,B,C,D,Ts,'zoh')
 %funcion de transferencia pulso Y(z)/U(z)
 [Nuz,Dez]=ss2tf(G,H,C,D)
 Gpz=tf(Nuz,Dez,Ts)
 
 %% Disño del Regulador METODO ALGORITMICO
 % Matrix se controlabilidad
 M=[H G*H G^2*H G^3*H]
 r=rank(M)
 % Rango = n=4 Completamente controlable
 %Polinomio de lazo abierto
 Paz=poly(G)
 a1=Paz(2);  a2=Paz(3); a3=Paz(4);  a4=Paz(5);
 W=[a3 a2 a1 1;a2 a1 1 0;a1 1 0 0;1 0 0 0]
 %Matrix de transformacion a su 1ra FCC
 T=M*W
 
 %Polos de lazo cerrado en tiempo continuo
 S=[-1+2j, -1-2j, -2, -4]
 %Polos de lazo cerrado en tiempo discreto
 Z=exp(S*Ts)
 %Polinomio de lazo cerrado
 Pcz=poly(Z)
 al1=Pcz(2);  al2=Pcz(3); al3=Pcz(4);  al4=Pcz(5);
 % Matriz de ganacia o ponderacion de estados
 K=[al4-a4  al3-a3  al2-a2  al1-a1]*inv(T)
 
 % Simulacion de sistema Regulador
 Gc=G-H*K
 Hc=[0;0;0;0]
 Cc=C-D*K
 Dc=0
 %Condiciones iniciales
 x0=[-0.2 0.3 0 0]'
 % Numero de muestras
 Nm=5/Ts+1
 k=0:Nm-1;
 [Y0, X0]=dinitial(Gc,Hc,Cc,Dc,x0,Nm);
 figure
 stairs(k*Ts,Y0)
 figure
 stairs(k*Ts,X0)
 legend('x1(k)','x2(k)','x3(k)','x4(k)')
 % Ley de Control U0
 U0=-K*X0';
 figure
 stairs(k*Ts,U0)
 
  
 %% Sistema de Seguimiento ante un escalon con ganancia en prealimentacion
 % Aplicando el metrodo algoritmico
 % Mp<=7%  ts=4 seg
 
  % Matrix se controlabilidad
 M=[H G*H G^2*H G^3*H]
 r=rank(M)
 % Rango = n=4 Completamente controlable
 %Polinomio de lazo abierto
 Paz=poly(G)
 a1=Paz(2);  a2=Paz(3); a3=Paz(4);  a4=Paz(5);
 W=[a3 a2 a1 1;a2 a1 1 0;a1 1 0 0;1 0 0 0]
 %Matrix de transformacion a su 1ra FCC
 T=M*W
 
 %Polos de lazo cerrado en tiempo continuo
 % a base de prueba y error
 S=[-1.2+1.1j, -1.2-1.1j, -2, -4]
 %Polos de lazo cerrado en tiempo discreto
 Z=exp(S*Ts)
 %Polinomio de lazo cerrado
 Pcz=poly(Z)
 al1=Pcz(2);  al2=Pcz(3); al3=Pcz(4);  al4=Pcz(5);
 % Matriz de ganacia o ponderacion de estados
 K=[al4-a4  al3-a3  al2-a2  al1-a1]*inv(T)
 % D=0 No se considera en la formula
 k0=inv(C*inv(eye(4)-G+H*K)*H)
 %Matrices de lazo cerrado
 Gc=G-H*K
 Hc=H*k0
 Cc=C-D*K
 Dc=D*k0
 Nm=7/Ts+1
 k=0:Nm-1;
 %Respuesta al escalon unitario SOLUCIO PARTICULAR C.I.=0
 [Yr, Xr]=dstep(Gc,Hc,Cc,Dc,1,Nm);
 figure
 stairs(k*Ts,Yr)
 figure
 stairs(k*Ts,Xr)
 legend('x1(k)','x2(k)','x3(k)','x4(k)')
 %Ley de Control Ur magitud del escalon r=1
 r=ones(Nm,1);
 Ur=-K*Xr'+k0*r';
 figure
 stairs(k*Ts,Ur)
 %Respuesta SOLUCIO HOMOGENEA CON C.I X0
 %Condiciones iniciales
 % Aplicacion del Regulador
 x0=[-0.2 0.3 0 0]'
 [Y0, X0]=dinitial(Gc,Hc*0,Cc,Dc*0,x0,Nm);
 U0=-K*X0';
 % SOLUCION GENERAL
 Y=X0+Xr;
 X=X0+Xr;
 U=U0+Ur;
 figure
 stairs(k*Ts,Y)
 figure
 stairs(k*Ts,X)
 legend('x1(k)','x2(k)','x3(k)','x4(k)')
 figure
 stairs(k*Ts,U) 
