%% Final Project Parts 1 and 2:  Reactive Ion Etching with MIMO Control

% Plant Model from [Power; Throttle] to [|F|; Vbias]
P1 = [tf([0.17 0.7],[1 15 26.7]); tf(0.28,[1 0.97])];
P2 = [tf(-0.17,[1 0.24]); tf([2.41 9.75],[1 4 0.7])];
P2.InputDelay = 0.5;
P = [P1 P2];
P_0 = freqresp(P,0);
num  = cond(P_0,2)
% Normalized System
DO = diag([30 350]);
DI = diag([1000 12.5]);
PN = inv(DO)*P*DI;
PN = ss(PN);

% Condition number of DC gain
PN_0 = freqresp(PN,0);
num  = cond(PN_0,2)

% Use second-order Pade approximation for input delay
PN = pade(PN,2);
PN.InputName = 'u';
PN.OutputName = 'y';

% State-space matrices and dimensions
[AP,BP,CP,DP] = ssdata(PN);
[nx,nu] = size(BP);
ny = size(CP,1);

%% Part 1(A): Linear Quadratic Regulator with Integrators

% Augment state equations so that you can do integral control
Aaug = [AP zeros(nx,ny);
    CP zeros(ny,ny)];
Baug = [BP;
    zeros(ny,nu)];
Caug = [CP zeros(ny,ny)];
%%
% LQR Weighting Matrices
Q1 = diag([2.5,2.5]);
Qi = diag([0.1,0.1]);
Q = CP'*Q1*CP;
R = diag([0.1,0.1]);
% LQ state feedback gain
K = lqr(Aaug,Baug,[Q zeros(nx,nu); zeros(nu,nx) Qi],R);
%%
% Closed loop state equations with state feedback and integrators.
sys = ss(Aaug-Baug*K,[zeros(nx,nu); -eye(2)],Caug,0);
sys1 = ss(Aaug-Baug*K,[zeros(nx,nu); -eye(2)],-K,0);
% Verify that closed-loop is stable (Check to verify no bugs in code)
isstable(sys);
% Time vector
Tf = 50;
Nt = 500;
t = linspace(0,Tf,Nt);

% Step responses
trans = tf(sys);
trans_1 = tf(sys1);

figure(1)
% set(findall(gcf,'type','line'),'linewidth',3);
subplot(221)
step(trans(1,1))
hold on
step(trans(2,1))
hold on
title('Step to [F], SF')
legend('[F]','Vbias')

subplot(222)
step(trans_1(1,1))
hold on
step(trans_1(2,1))
hold on
title('Step to [F], SF')
legend('RF Power','Throttle')

subplot(223)
step(trans(1,2))
hold on
step(trans(2,2))
hold on
title('Step to Vbias, SF')
legend('[F]','Vbias')

subplot(224)
step(trans_1(1,2))
hold on
step(trans_1(2,2))
hold on
title('Step to Vbias, SF')
legend('RF Power','Throttle')
%% Part 1(B): Linear Quadratic Regulator with Integrators

% Covariance matrices for loop transfer recovery observer
q = 2.5;
V =  q^2*BP*BP';
W = eye(2);

% Observer gain
% Note: We only need to estimate the plant states. We do not need the
% observer to construct an estimate of the integrator states.
G = eye(8);
L = lqe(AP,G,CP,V,W);
Laug = [L;
    eye(2)];
Cobs = ss([AP-BP*K(:,1:8)-L*CP -BP*K(:,9:10); zeros(2,10)],Laug,K,0);
% sys3 = ss([AP-L*CP-BP*K(:,1:8) -BP*K(:,9:10); zeros(2,8) zeros(2,2)],)
loop_sys3 = Cobs*PN;
sys3 = feedback(loop_sys3,eye(2));
sys4 = tf(Cobs)*feedback(eye(2),loop_sys3);
%%
% Verify that observer error is stable (Check to verify no bugs in code)
eig2 = eig(Aaug-Baug*K);
eig3 = eig(AP-L*CP);

% Construct controller: 
% This includes the observer, integrators, and feedback gains.
%   Inputs: [|F| Ref.; Vbias Ref; |F| measurement; Vbias measurement]
%   Outputs: [Power; Throttle]
A_obs = [AP -BP*K(:,1:8) -BP*K(:,9:10);
        L*CP AP-BP*K(:,1:8)-L*CP -BP*K(:,9:10);
        CP zeros(2,10)];
B_obs = [[zeros(16,2); eye(2)] [zeros(8,1); L(:,1); zeros(2,1)] [zeros(8,1); L(:,1); zeros(2,1)]];
C_obs= [CP zeros(size(CP)) zeros(2,2); zeros(2,8) -K];

% Form Closed-Loop
%   Inputs are [|F| Ref.; Vbias Ref; |F| noise; Vbias noise]
%   Outputs: [|F|; Vbias; Power; Throttle]
sys5 = -ss(A_obs,B_obs,C_obs,0);
sys6 = ss(A_obs,B_obs,[zeros(2,8) K],0);

% Verify that closed-loop eigenvalues are the union of the observer
% and state-feedback eigenvalues. This is a useful debugging step
% to verify that you have correctly formed the closed-loop.
eig1 = eig(A_obs);
isstable(sys5);

% trans = tf(sys3);
% trans_1 = tf(sys4);
%%
% Step responses with noise on |F|
noise=0.1*randn(501,1);
ref1 = [ones(size(noise)) zeros(size(noise)) noise zeros(size(noise))];
noise=0.1*randn(501,1);
ref2 = [zeros(size(noise)) ones(size(noise)) zeros(size(noise)) noise];
[y1 t]=lsim(sys5,ref1,[0:0.1:50]);

[y2 t]=lsim(sys5,ref2,[0:0.1:50]);

figure(1)
subplot(221)
plot(t,y1(:,1),t,y1(:,2))
hold on
title('Step to [F], MIMO')
legend('[F]','Vbias')
xlim([0 50])

subplot(222)
plot(t,y1(:,3),t,y1(:,4))
hold on
title('Step to [F], MIMO')
legend('RF Power','Throttle')
xlim([0 50])

subplot(223)
plot(t,y2(:,1),t,y2(:,2))
hold on
title('Step to [F], MIMO')
legend('[F]','Vbias')
xlim([0 50])

subplot(224)
plot(t,y2(:,3),t,y2(:,4))
hold on
title('Step to [F], MIMO')
legend('RF Power','Throttle')
xlim([0 50])
%%
% Bode magnitude from |F| noise to [|F|; Vbias; Power; Throttle]
figure(5)
bodemag(sys5(1,3),sys5(2,3),sys5(3,3),sys5(4,3))
legend('F ','Vbias','RF Power','Throttle')
% Sigma magnitude from [|F| Ref.; Vbias Ref] to [Power; Throttle]
figure(6)
sigma(sys5(3,1),sys5(3,2),sys5(4,1),sys5(4,2))
title('Sigma magnitude from [|F| Ref.; Vbias Ref] to [Power; Throttle]')
legend('|F| Ref. to Power',' Vbias Ref to Power','|F| Ref. to Throttle',' Vbias Ref to Throttle')
% Sigma magnitude from [|F| Noise; Vbias Noise] to [Power; Throttle]
figure(7)
sigma(sys5(3,3),sys5(3,4),sys5(4,3),sys5(4,4))
title('Sigma magnitude from [|F| Noise; Vbias Noise] to [Power; Throttle]')
legend('|F| Noise to Power',' Vbias Noise to Power','|F| Noise to Throttle',' Vbias Noise to Throttle')

%% Part 1(D): Stability Margins and Comparison With State Feedback

% Loop-at-a-time margins at the plant input
LI = (Cobs * PN);

AM = allmargin (LI)
AM(1)%For the first channel
AM(2)%For the second channel

Lsf = ss(Aaug,Baug,K,[])
% Multi-Loop Disk margins at the plant input
[DMI_sf,MMI_sf] = diskmargin(Lsf);

[DMI,MMI] = diskmargin(Cobs*PN);  % T-based disk margin
MMI

% Unstructured (fully-coupled) stability margin (USM) at the plant input
TI = feedback(LI,eye(size(LI)));
[np,wp] = hinfnorm(TI);
StabMarg1 = 1/np % unstructured stability margin

% Input loop transfer function: Compare Lsf to LI
Lsf = ss(Aaug,Baug,K,[])
Tsf = feedback(Lsf,eye(size(Lsf)));
TI = feedback(LI,eye(size(LI)));


figure (5) % for Loop Gain
[sv,wout]=sigma(LI,{1e-2, 1e2});
loglog(wout,sv(1,:),wout,sv(2,:))
grid on
hold on
[sv,wout]=sigma(Lsf,{1e-2, 1e2});
loglog(wout,sv(1,:),wout,sv(2,:))
hold off
title('Open loop singular values LsF(s) = K*inv(sI-A)*B and Lobs(s) = C(s)*P(s)')
legend('sig_{max}(Lobs)','sig_{min}(Lobs)','sig_{max}(LsF)','sig_{min}(LsF)')
xlabel('w rad/sec')
ylabel('Magnitude')

figure(6) % for Complimentary sensitivity
[sv,wout]=sigma(TI,{1e-2, 1e2});
loglog(wout,sv(1,:),wout,sv(2,:))
grid on
hold on
[sv,wout]=sigma(Tsf,{1e-2, 1e2});
loglog(wout,sv(1,:),wout,sv(2,:))
hold off
title('singular values of complementary sensitivity function')
legend('sig_{max}(Tobs)','sig_{min}(Tobs)','sig_{max}(TsF)','sig_{min}(TsF)')
xlabel('w rad/sec')
ylabel('Magnitude')
Ssf = feedback(eye(size(Lsf)),Lsf);
SI = feedback(eye(size(LI)),LI);

figure (7) % for Sensitivity
[sv,wout]=sigma(SI,{1e-2, 1e2});
loglog(wout,sv(1,:),wout,sv(2,:))
grid on
hold on
[sv,wout]=sigma(Ssf,{1e-2, 1e2});
loglog(wout,sv(1,:),wout,sv(2,:))
hold off
title('singular values of sensitivity function')
legend('sig_{max}(Sobs)','sig_{min}(Sobs)','sig_{max}(SsF)','sig_{min}(SsF)')
xlabel('w rad/sec')
ylabel('Magnitude')
% figure (11)
% y1= (sigma(db2mag(TI,{1e-2, 1e2})))
% y2= (sigma(db2mag(Tsf,{1e-2, 1e2})))

% plot({1e-2, 1e2},y1,{1e-2, 1e2},y2)

% Input sensitivity: Compare Ssf to SI

% Input complementary sensitivity: Compare Tsf to TI


%% Part 2(B): Equivalent Controller

% Equivalent controller
%   Ceq = inv[ I+K1 inv(sI-A) B] (KI/s)
sys = ss(AP,BP,K(:,1:8),0);
tr_ = tf(1,[1 0]);
% Ceq = tf(inv( eye(2) + K(:,1:8)*inv(s*eye(8) - AP)*BP)*K(:,9:10)/s);
Ceq = inv(eye(2) + sys)*K(:,9:10)*tr_;
figure(8)
bode(Ceq(1,1))
hold on
bode(Ceq(1,2))
hold on
bode(Ceq(2,1))
hold on
bode(Ceq(2,2))
hold on
legend('Ceq(1,1)','Ceq(1,2)','Ceq(2,1)','Ceq(2,2)')
title('equivalent compensator')

%% Part 2(C): Decentralized Approximation of Equivalent Controller
s = tf('s');
Chateq1 = 0.775*((s+13)*(s+0.25)*(s+0.97)*(s+2)/((s+31)*(s+4.1)*(s+1.2)*s*(s+0.2)));
% Chateq2 = 0.632*((s+13)*(s+2)*(s+0.96)*(s+0.97)*(s+0.18)/((s+31)*(s)*(s+4.1)*(s+0.9)*(s+0.21)*(s+1.17)));
Chateq2 = 0.775*((s+0.24)*(s+0.19)/((s+0.9)*(s)*(s+0.21)))
% Chateq2 = tf([1],[1 0])*tf(0.38*[1 .203],0.203*[1 0.38])*tf([0.1275]);
Chateq = [Chateq1 Chateq1; -Chateq2 Chateq2];
figure(9)
bode(Ceq(1,1),10^-2:10^-2:10^3)
hold on
bode(Chateq(1,1),10^-2:10^-3:10^3)
hold on
bode(Ceq(2,2),10^-2:10^-2:10^3)
bode(Chateq(2,2),10^-2:10^-2:10^3)
hold on
legend('Ceq(1,1)','Cd1','Ceq(2,2)','Cd2')

figure(10)
step(feedback(PN*Ceq,eye(2)))
hold on
step(feedback(PN*Chateq,eye(2)))
hold on
legend('with Ceq','with approx to Ceq')
title('step response with Ceq and approximation')
%% Part 2: Comment on Plant Transformation
M = [1 1; -1 1]/sqrt(2);
MP = M*PN;

figure(12)
subplot(2,1,1)
bodemag(PN(1,1),'b',PN(1,2),'r--',PN(2,1),'m-.',PN(2,2),'g-.',{1e-2,1e2});
legend('PN(1,1)','PN(1,2)','PN(2,1)','PN(2,2)','Location','Southwest');
grid on;
if exist('garyfyFigure','file'), garyfyFigure, end

subplot(2,1,2)
bodemag(MP(1,1),'b',MP(1,2),'r--',MP(2,1),'m-.',MP(2,2),'g-.',{1e-2,1e2});
legend('MP(1,1)','MP(1,2)','MP(2,1)','MP(2,2)','Location','Southwest');
grid on;
if exist('garyfyFigure','file'), garyfyFigure, end

