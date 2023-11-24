%%%%%%%%%%%% deterministic signals%%%%%%%%%%%%
%%%%% periodic waveform %%%%%%%%%%%%%%%%%%%%%%
%%%% https://fr.mathworks.com/help/signal/ug/signal-generation-and-visualization.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 10000;
t = 0:1/fs:1.5;
x1 = sawtooth(2*pi*50*t);
x2 = square(2*pi*50*t);

figure
nexttile
plot(t,x1)
axis([0 0.2 -1.2 1.2])
xlabel("Time (sec)")
ylabel("Amplitude") 
title("Sawtooth Periodic Wave")

nexttile
plot(t,x2)
axis([0 0.2 -1.2 1.2])
xlabel("Time (sec)")
ylabel("Amplitude")
title("Square Periodic Wave")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% aperiodic waveform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 10000;
t = -1:1/fs:1;
x1 = tripuls(t,20e-3);
x2 = rectpuls(t,20e-3);

figure
nexttile
plot(t,x1)
axis([-0.1 0.1 -0.2 1.2])
xlabel("Time (sec)")
ylabel("Amplitude")
title("Triangular Aperiodic Pulse")

nexttile
plot(t,x2)
axis([-0.1 0.1 -0.2 1.2])
xlabel("Time (sec)")
ylabel("Amplitude")
title("Rectangular Aperiodic Pulse")


tc = gauspuls("cutoff",50e3,0.6,[],-40); 
t1 = -tc : 1e-6 : tc; 
y1 = gauspuls(t1,50e3,0.6);

nexttile
plot(t1,y1)
axis([-0.1 0.1 -0.2 1.2])
xlabel("Time (sec)")
ylabel("Amplitude")
title("Gaussian  Aperiodic Pulse")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Swept-Frequency Waveforms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tlin = 0:0.001:2;
ylin = chirp(tlin,100,1,250);

tq = -2:0.001:2;
yq = chirp(tq,100,1,200,"quadratic");

nexttile
plot(tlin,ylin)
axis([0 2 -0.2 1.2])
xlabel("Time (sec)")
ylabel("Amplitude")
title("Chirp")

nexttile
plot(tq,yq)
axis([0 2 -0.2 1.2])
xlabel("Time (sec)")
ylabel("Amplitude")
title("Quadratic Chirp")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% pulse train %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 100e9;
D = [2.5 10 17.5]' * 1e-9;
t = 0 : 1/fs : 2500/fs;
w = 1e-9;
yp = pulstran(t,D,@rectpuls,w);
T = 0 : 1/50e3 : 10e-3;
D = [0 : 1/1e3 : 10e-3 ; 0.8.^(0:10)]';
Y = pulstran(T,D,@gauspuls,10e3,.5);

figure
nexttile
plot(t*1e9,yp)
axis([0 25 -0.2 1.2])
xlabel("Time (ns)")
ylabel("Amplitude")
title("Rectangular Train")

nexttile
plot(T*1e3,Y)
xlabel("Time (ms)")
ylabel("Amplitude")
title("Gaussian Pulse Train")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  DFA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=8;n=1;
[x,t,Fs]=sinus(f,n);
pts = 5:5:(size(x,2)-1);
[A,F,Y,Y_hat] = DFA_function(x,pts,1);
plot_fun = @(xp,A,ord) polyval(A,log(xp));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(2,4,1);
title('Regression');
plot(Y(3,1:100),'r'),
hold on,
plot(Y_hat(3,1:100),'b'),
hold off
xlabel('fenêtre: 45');
subplot(2,4,2);
plot(Y(5,:),'r'),
hold on,
plot(Y_hat(5,:),'b'),
hold off
subplot(2,4,3);
plot(Y(15,:),'r'),
hold on,
plot(Y_hat(15,:),'b'),
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
title('DFA sinusoïde sur une période');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=8;n=50;
[x,t,Fs]=sinus(f,n);
pts = 5:20:(size(x,2)-1);
[A,F,Y,Y_hat] = DFA_function(x,pts,1);
figure(2)
subplot(2,4,1);
title('DFA sinusoïde sur 50 période')
plot(Y(3,:),'r'),
hold on,
plot(Y_hat(3,:),'b'),
hold off
xlabel('fenêtre: 45');
subplot(2,4,2);
plot(Y(20,:),'r'),
hold on,
plot(Y_hat(20,:),'b'),
hold off
subplot(2,4,3);
plot(Y(50,:),'r'),
hold on,
plot(Y_hat(50,:),'b'),
hold off
subplot(2,4,4);
scatter(log(pts),log(F));
title('DFA sinusoïde 50 périodes');
hold  on
plot(log(pts),plot_fun(pts,A,'--'));
xlabel(A(1));
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=8;n=50;
[x,t,Fs]=sinus(f,n);
pts = 5:20:(size(x,2)-1);
[A_1,F_1,Y_1,Y_hat_1] = DFA_function(x,pts(1:5),1);
[A_2,F_2,Y_2,Y_hat_2] = DFA_function(x,pts(6:end),1);

figure(3)
subplot(2,4,1);
plot(Y_1(3,1:200),'r'),
hold on,
plot(Y_hat_1(3,1:200),'b'),
hold off
subplot(2,4,2);
plot(Y_1(5,1:500),'r'),
hold on,
plot(Y_hat_1(5,1:500),'b'),
hold off
subplot(2,4,3);
plot(Y_2(50,1:5000),'r'),
hold on,
plot(Y_hat_2(50,1:5000),'b'),
hold off
subplot(2,4,4);
scatter(log(pts(1:5)),log(F_1(1:5)));
title('DFA');
hold  on
plot(log(pts(1:5)),plot_fun(pts(1:5),A_1,'--'));
xlabel(A_1(1));
scatter(log(pts(6:end)),log(F_2));
title('DFA sinusoïde 50 périodes');
plot(log(pts(6:end)),plot_fun(pts(6:end),A_2,'--'));
xlabel(A_1(1),A_2(1));
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear all;
    f2=100,
    f1=8;n=10;
    [x2,t2,Fs]=sinus2f(f1,f2,n);
    

    pts_1=2:2:(size(x2,2)-1);
    pts_2=2:8:(size(x2,2)-1);
    pts_3=2:32:(size(x2,2)-1);
    pts_4=2:128:(size(x2,2)-1);
    pts_5=2:400:(size(x2,2)-1);
    e = 1:10:1000;
    plot_fun = @(xp,A,ord) polyval(A,log(xp));

    [Asin21,Fsin21,Y_21,Y_hat_21] = DFA_function(x2,pts_1,1);
    [Asin22,Fsin22,Y_22,Y_hat_22] = DFA_function(x2,pts_2,1);
    [Asin23,Fsin23,Y_23,Y_hat_23] = DFA_function(x2,pts_3,1);
    [Asin24,Fsin24,Y_24,Y_hat_24] = DFA_function(x2,pts_4,1);
    [Asin25,Fsin25,Y_25,Y_hat_25] = DFA_function(x2,pts_5,1);


    figure(4)
    subplot(2,4,1);
    scatter(log(pts_1),log(Fsin21))
    title('DFA, pas 2')
    hold on
    plot(log(e),plot_fun(e,Asin21),'--')
    xlabel(Asin21(1));
    subplot(2,4,2); 
    scatter(log(pts_2),log(Fsin22))
    title('DFA, pas 8');
    hold on
    plot(log(e),plot_fun(e,Asin22),'--')
    xlabel(Asin22(1));
    subplot(2,4,3);
    scatter(log(pts_3),log(Fsin23))
    title('DFA, pas 32');
    hold on
    plot(log(e),plot_fun(e,Asin23),'--')
    xlabel(Asin23(1));
    subplot(2,4,4);
    scatter(log(pts_5),log(Fsin25))
    title('DFA sin f1,f2, pas 128');
    hold on
    plot(log(e),plot_fun(e,Asin25),'--')
    xlabel(Asin25(1));
    hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 10000;
t = 0:1/fs:1.5;
xs = sawtooth(2*pi*50*t);
xq = square(2*pi*50*t);
pts = 5:20:(size(x1,2)-1);
[A,F,Y,Y_hat] = DFA_function(x1,pts,1);
plot_fun = @(xp,A,ord) polyval(A,log(xp));
figure(1)
subplot(2,4,1);
title('Regression Square');
plot(Y(5,:),'r'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(5,:),'b'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold off
xlabel('time');
subplot(2,4,2);
plot(Y(10,:),'r'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(10,:),'b'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(100,:),'r'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(100,:),'b'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
title('DFA sawtooth');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off

[A,F,Y,Y_hat] = DFA_function(xq,pts,1);
figure(2)
subplot(2,4,1);
title('Regression');
plot(Y(5,:),'r'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(5,:),'b'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold off
xlabel('time');
subplot(2,4,2);
plot(Y(10,:),'r'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(10,:),'b'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(100,:),'r'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(100,:),'b'),
axis([0 500 min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
title('DFA square wave');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% aperiodic waveform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
fs = 10000;
t = -1:1/fs:1;
x1 = tripuls(t,20e-3);
x2 = rectpuls(t,20e-3);

pts = 5:20:(size(x1,2)-1);
[A,F,Y,Y_hat] = DFA_function(x1,pts,1);
plot_fun = @(xp,A,ord) polyval(A,log(xp));


figure(1)
subplot(2,4,1);

plot(Y(30,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(30,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,2);
plot(Y(400,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(400,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(900,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(900,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))

[A,F,Y,Y_hat] = DFA_function(x2,pts,1);
plot_fun = @(xp,A,ord) polyval(A,log(xp));

figure(2)
subplot(2,4,1);
plot(Y(30,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(30,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,2);
plot(Y(400,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(400,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(900,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(900,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
title('DFA aperiodic tripuls');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off
title('DFA aperiodic rectpuls');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off


tc = gauspuls("cutoff",50e3,0.6,[],-40); 
t = -tc : 1e-6 : tc; 
y1 = gauspuls(t,50e3,0.6);

[A,F,Y,Y_hat] = DFA_function(y1,pts,1);
figure(2)
subplot(2,4,1);
plot(Y(1,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(1,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,2);
plot(Y(3,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(3,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(900,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(900,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
title('DFA aperiodic tripuls');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off
title('DFA gaussian pulse');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Swept-Frequency Waveforms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tlin = 0:0.001:2;
ylin = chirp(tlin,100,1,250);
pts = 5:20:(size(ylin,2)-1);
tq = -2:0.001:2;
yq = chirp(tq,100,1,200,"quadratic");
t=tlin;
[A,F,Y,Y_hat] = DFA_function(ylin,pts,1);
plot_fun = @(xp,A,ord) polyval(A,log(xp));
figure(2)
subplot(2,4,1);
plot(Y(3,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(3,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,2);
plot(Y(20,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(20,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(90,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(90,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
title('DFA aperiodic tripuls');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off
title('DFA sweep chirp');
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
hold off

t=tq;
[A,F,Y,Y_hat] = DFA_function(yq,pts,1);
figure(2)
subplot(2,4,1);
plot(Y(3,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(3,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,2);
plot(Y(20,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(20,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(90,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(90,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
title('DFA sweep chirp quadratic');
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% pulse train %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 100e9;
D = [2.5 10 17.5]' * 1e-9;
t = 0 : 1/fs : 2500/fs;
w = 1e-9;
yp = pulstran(t,D,@rectpuls,w);
T = 0 : 1/50e3 : 10e-3;
D = [0 : 1/1e3 : 10e-3 ; 0.8.^(0:10)]';
yg = pulstran(T,D,@gauspuls,10e3,.5);


pts = 5:20:(size(yp,2)-1);
[A,F,Y,Y_hat] = DFA_function(yp,pts,1);
plot_fun = @(xp,A,ord) polyval(A,log(xp));
figure(6)
subplot(2,4,1);
plot(Y(3,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(3,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,2);
plot(Y(20,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(20,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(90,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(90,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
title('DFA pulse train');
hold off


pts = 5:20:(size(yg,2)-1);
[A,F,Y,Y_hat] = DFA_function(yg,pts,1);
plot_fun = @(xp,A,ord) polyval(A,log(xp));
t=T;
figure(6)
subplot(2,4,1);
plot(Y(3,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(3,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,2);
plot(Y(10,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(10,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(20,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(20,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
title('DFA Gaussian pulse train');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Noise         %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pts = 50:20:2000;

for i=1:50
    wn = dsp.ColoredNoise('white',4000);
    pn = dsp.ColoredNoise('Pink',4000);
    bn = dsp.ColoredNoise('brown',4000);
    Wn = wn();
    Pn = pn();
    Bn = bn();
    [Aw,Fw] = DFA_function(Wn,pts,1);
    [Ap,Fp] = DFA_function(Pn,pts,1);
    [Ab,Fb] = DFA_function(Bn,pts,1);
    Alpha_w(i)=Aw(1);
    Alpha_p(i)=Ap(1);
    Alpha_b(i)=Ab(1);
end;
figure(1),
plot(Alpha_w,'k');
hold on
plot(Alpha_p,'r');
title('Alpha 50 itérations')
hold on
plot(Alpha_b,'b');
legend('wite noise','pink noise','Brown noise')
hold off 



bn = dsp.ColoredNoise('brown',20000);
Bn=bn();
pts = 5:20:(size(Bn',2)-1);
[A,F,Y,Y_hat] = DFA_function(Bn,pts,1);
plot_fun = @(xp,A,ord) polyval(A,log(xp));
t=Bn';
figure(6)
subplot(2,4,1);
plot(Y(3,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(3,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,2);
plot(Y(100,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(100,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,3);
plot(Y(200,:),'r'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold on,
plot(Y_hat(200,:),'b'),
axis([0 size(t,2) min(Y(1,:)) max(Y(1,:))])
hold off
subplot(2,4,4);
scatter(log(pts),log(F))
hold on
plot(log(pts),plot_fun(pts,A),'--')
xlabel(A(1));
title('DFA Brown Noise');
hold off

%%
% 
% <html>
% <table border=1><tr><td>one</td><td>two</td></tr></table>
% </html>
% 



