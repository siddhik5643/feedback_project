%% Final Project Part 3:  Reactive Ion Etching with MIMO Control

%% Part 3(A) -- DC Analysis With Oxygen Sensor

% Model
% Inputs: [Power; Throttle; %O2]
% Outputs: [|F|; Vbias; Pressure]
P1 = [zpk(-0.067,[-0.095 -19.69],0.49); ...
    zpk(-0.27,[-0.19 -62.42],12.23); ...
    zpk(0.006,[-0.19 -2.33],-0.011)];

P2 = [zpk(0.73,[-0.11; -39.76],4.85); tf(1.65,[1 0.16]); ...
      zpk([],[-0.18; -3],-0.97)];
P2.InputDelay = 0.42;

P3 = [tf(0.33,[1 0.17]); tf(0.25,[1 0.41]); tf(0.024,[1 0.4])];
P3.InputDelay = 0.77;

Pox = [P1 P2 P3];

% Use second-order Pade for plant
Pox = pade(Pox, 2);

% Condition number of DC gain for not scaled plant 
Pox_0 = freqresp(Pox,0);
num  = cond(Pox_0,2)

% Input and output scalings based on equilibrium values
DO = diag([16.52 340 17.83]);
DI = diag([1000 12.5 5]);

% Normalize Plant
PoxN = inv(DO)*Pox*DI;
PoxN = ss(PoxN);
% Condition number of DC gain
PoxN_0 = freqresp(PoxN,0);
num  = cond(PoxN_0,2)

%% Part 3(B) -- DC Analysis With Oxygen Sensor

% State-space data for scaled plant
[As,Bs,Cs] = ssdata(PoxN);
[nx,nu] = size(Bs);
ny = size(Cs,1);

% Weighting matrices (Q,R,V,W)
% Assume Q of the form Q = blkdiag(alpha*Cs'*Cs, Qw)

% Q1 = diag([1,1,1]);
% Qi = diag([0.1,0.1,0.1]);
% Q = blkdiag(Cs'*Q1*Cs,Qi);
% R = diag([1,1,1]);

Q1 = diag([2.5,2.5,2.5]);
Qi = diag([0.2,0.2,0.2]);
Q = blkdiag(Cs'*Q1*Cs,Qi);
R = diag([0.1,0.1,0.2]);

% Augmented Plant with integrators
Aaug = [As zeros(nx,ny);
    Cs zeros(ny,ny)];
Baug = [Bs;
    zeros(ny,nu)];
Caug = [Cs zeros(ny,ny)];

% Compute state feedback and observer gains
K = lqr(Aaug,Baug,Q,R);
q = 3;
V =  q^2*Bs*Bs';
W = eye(3);
G = eye(26);
L = lqe(As,G,Cs,V,W);

sys = ss(Aaug-Baug*K,[zeros(nx,nu); -eye(nu)],Caug,0);
sys1 = ss(Aaug-Baug*K,[zeros(nx,nu); -eye(nu)],-K,0);
% %%
% % Time vector
% Tf = 50;
% Nt = 500;
% t = linspace(0,Tf,Nt);
% 
% % Step responses
% trans = tf(sys);
% trans_1 = tf(sys1);
% 
% figure(1)
% 
% subplot(2,3,1)
% step(trans(1,1))
% hold on
% step(trans(2,1))
% hold on 
% step(trans(3,1))
% title('Step to [F], SF')
% grid on
% 
% subplot(2,3,2)
% step(trans(1,2))
% hold on
% step(trans(2,2))
% hold on 
% step(trans(3,2))
% title('Step to [F], SF')
% grid on
% 
% subplot(2,3,3)
% step(trans(1,3))
% hold on
% step(trans(2,3))
% hold on 
% step(trans(3,3))
% title('Step to Vbias, SF')
% 
% grid on
% 
% subplot(2,3,4)
% step(trans_1(1,1))
% hold on
% step(trans_1(2,1))
% hold on 
% step(trans_1(3,1))
% title('Step to Vbias, SF')
% 
% grid on 
% 
% subplot(2,3,5)
% step(trans_1(1,2))
% hold on
% step(trans_1(2,2))
% hold on 
% step(trans_1(3,2))
% 
% 
% subplot(2,3,6)
% step(trans_1(1,3))
% hold on
% step(trans_1(2,3))
% hold on 
% step(trans_1(3,3))
%%
% Construct controller: 
% This includes the observer, integrators, and feedback gains.
%   Inputs: [|F|Ref; Vbias Ref; Press Ref; |F| Meas; Vbias Meas; Press Meas]
%   Outputs: [Power; Throttle; %O2]
A_obs = [As -Bs*K(:,1:26) -Bs*K(:,27:29);
        L*Cs As-Bs*K(:,1:26)-L*Cs -Bs*K(:,27:29);
        Cs zeros(3,29)];
B_obs = [[zeros(52,3); eye(3)] [zeros(26,1); L(:,1); zeros(3,1)] [zeros(26,1); L(:,1); zeros(3,1)] [zeros(26,1); L(:,1); zeros(3,1)]];
C_obs= [Cs zeros(size(Cs)) zeros(3,3); zeros(3,26) -K];
sys5 = -ss(A_obs,B_obs,C_obs,0);
sys6 = ss(A_obs,B_obs,[zeros(3,26) K],0);

% Form Closed-Loop
%   Inputs: [|F|Ref; Vbias Ref; Press Ref; |F| Noise; Vbias Noise; Press Noise]
%   Outputs: [|F|; Vbias; Press; Power; Throttle; %O2]
trans = tf(sys5);
trans_1 = tf(sys6);


% Verify that closed-loop eigenvalues are the union of the observer
% and state-feedback eigenvalues. This is a useful debugging step
% to verify that you have correctly formed the closed-loop.
eig2 = eig(Aaug-Baug*K);
eig3 = eig(AP-L*CP);
eig1 = eig(A_obs);
isstable(sys5);

%%
% Time vector
Tf = 40;
Nt = 400;
t = linspace(0,Tf,Nt);

% Step Responses (Without Noise)
noise=0.1*randn(401,1);
ref1 = [ones(size(noise)) zeros(size(noise)) zeros(size(noise)) zeros(size(noise)) zeros(size(noise)) zeros(size(noise))];
noise=0.1*randn(401,1);
ref2 = [zeros(size(noise)) ones(size(noise)) zeros(size(noise)) zeros(size(noise)) zeros(size(noise)) zeros(size(noise))];
ref3 = [zeros(size(noise)) zeros(size(noise)) ones(size(noise)) zeros(size(noise)) zeros(size(noise)) zeros(size(noise))];
[y1 t]=lsim(sys5,ref1,[0:0.1:40]);

[y2 t]=lsim(sys5,ref2,[0:0.1:40]);

[y3 t]=lsim(sys5,ref3,[0:0.1:40]);

figure

subplot(231)
plot(t,y1(:,1),t,y1(:,2),t,y1(:,3))
title('Step to [F]')
grid on

subplot(234)
plot(t,y1(:,4),t,y1(:,5),t,y1(:,6))
title('Step to [F]')
grid on

subplot(232)
plot(t,y2(:,1),t,y2(:,2),t,y2(:,3))
title('Step to Vbias')

grid on

subplot(235)
plot(t,y2(:,4),t,y2(:,5),t,y2(:,6))

title('Step to Vbias')

grid on

subplot(233)
plot(t,y3(:,1),t,y3(:,2),t,y3(:,3))
title('Step to Pressure')
legend('[F]','Vbias','Pressue')
grid on

subplot(236)
plot(t,y3(:,4),t,y3(:,5),t,y3(:,6))
title('Step to Pressure')
legend('RF Power','Throttle','%O2')
grid on 
