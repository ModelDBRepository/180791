function twofm()
% function twofm()
% A code used to simulate the two-factor model of OD plasticity 
% March 2015 by Taro Toyoizumi
%
% Reference:
% Toyoizumi T, Kaneko M, Stryker MP, Miller KD (2014)
% Modeling the dynamic interaction of hebbian and homeostatic plasticity.
% Neuron 84:497-510

%%% Parameters %%%
rcontra=.62; %ratio of contra-eye neurons
P.N.total=500; %total #neurons 
P.N.contra=round(rcontra*P.N.total); %#contra-eye neurons 
P.N.ipsi=P.N.total-P.N.contra; %#ipsi-eye neurons 
P.Input.mu=1.; %baseline firing rate
P.Input.BE=.5; %between-eye correlation
P.Input.MD=.5; %MD factor
P.Input.L=.2; %input correlation width
P.th.corr=.6; %constant threshold
P.y0=1.; %homeostatic setpoint
P.rhom=.7; %minimum value of Hebbian factor
P.rhoM=1.; %maximum value of Hebbian factor
P.tau.rho=.2; %Hebbian time-constant (day) 
P.tau.Hh=4.; %homeostatic hidden averaging time (day) 
P.T.wu=50; %warm-up (day)  
P.T.baseline=1; %warm-up (day)  
P.T.MD=7; %MD period (day)  
P.T.recovery=3; %recovery period (day)  
P.Hfun = @(x) max(1,x); 
P.Hhfun = @(x) (1+tanh(x-1)).*(x>1.05);
flag_trk=1;
flag_tnf=1;
flag_nmda=1;

%%% Simulation %%%
stages=fieldnames(P.T);
[Mu,Q,S]=InputStat(P.N,P.Input,0);

Init.time=-P.T.wu;
Init.rho=P.rhom+(P.rhoM-P.rhom)*[abs((1:P.N.contra)'-.5*P.N.contra)<.25*P.N.contra; abs((1:P.N.ipsi)'-.5*P.N.ipsi)<.25*P.N.ipsi];
Init.Hh=0;
Init.H=1;
    
Flags(1)=struct('MD',0,'Trk',1,'TNF',1,'NMDA',1);
Flags(2)=struct('MD',0,'Trk',1,'TNF',flag_tnf,'NMDA',1);
Flags(3)=struct('MD',1,'Trk',flag_trk,'TNF',flag_tnf,'NMDA',flag_nmda);
Flags(4)=struct('MD',0,'Trk',flag_trk,'TNF',flag_tnf,'NMDA',1);
    
Out0=PlastSim(Init,S,Flags(1),P,stages{1}); %warm-up
Hist=FinalValues(Out0);
for i=2:length(stages)
    disp(stages{i});
    Out=PlastSim(FinalValues(Hist),S,Flags(i),P,stages{i});
    Hist=AppendStruct(Hist,Out);
end

%%% Plotting %%%
FS=15;
Hist.time=Hist.time-P.T.baseline;
[Time,Ind]=meshgrid(Hist.time,1:P.N.total);
figure(1); clf;

% panel 1
subplot(2,3,1); 
pcolor(Time,Ind,Hist.w); cb=colorbar; shading flat;
set(gca,'fontsize',FS,'xtick',[0,P.T.MD,P.T.MD+P.T.recovery],'ytick',[1,round(.5*P.N.contra),P.N.contra,P.N.contra+round(.5*P.N.ipsi),P.N.total],'yticklabel',{'-','C','-','I','-'});
title(cb,'w'); xlabel('time (day) '); ylabel('index');

% panel 2
subplot(2,3,2); 
pcolor(Time,Ind,Hist.rho); shading flat; cb=colorbar; 
set(gca,'fontsize',FS,'xtick',[0,P.T.MD,P.T.MD+P.T.recovery],'ytick',[1,round(.5*P.N.contra),P.N.contra,P.N.contra+round(.5*P.N.ipsi),P.N.total],'yticklabel',{'-','C','-','I','-'});
title(cb,'r'); xlabel('time (day) '); ylabel('index');

% panel 3
subplot(2,3,3); 
plot(Hist.time,Hist.ODI,'k-','linewidth',2);
set(gca,'fontsize',FS,'xtick',[0,3,P.T.MD,P.T.MD+P.T.recovery]); xlim([-P.T.baseline,P.T.MD+P.T.recovery]); grid on;
xlabel('time (day) '); ylabel('ODI'); ylim([0 .4]);

% panel 4
subplot(2,3,4); 
indC=round(.5*P.N.contra);
indI=P.N.contra+round(.5*P.N.ipsi);
plot(Hist.time,Hist.resp(1,:)./Hist.resp(1,1),'r-',Hist.time,Hist.resp(2,:)./Hist.resp(2,1),'b-','linewidth',2); grid on;
set(gca,'fontsize',FS,'xtick',[0,3,P.T.MD,P.T.MD+P.T.recovery]); xlim([-P.T.baseline,P.T.MD+P.T.recovery]); 
ylim([.6,1.4]);
xlabel('time (day) '); ylabel('response ');

% panel 5
subplot(2,3,5); 
plot(Hist.time,Hist.Hh,'g-',Hist.time,Hist.H,'k-','linewidth',2); xlim([-P.T.baseline,P.T.MD+P.T.recovery]); 
set(gca,'fontsize',FS,'xtick',[0,3,P.T.MD,P.T.MD+P.T.recovery]); grid on;
xlabel('time (day) '); ylabel('H, h'); ylim([-.02,1.5]); 
    
% panel 6
subplot(2,3,6);
plot(Hist.time,Hist.y,'k-','linewidth',2);
set(gca,'fontsize',FS,'xtick',[0,3,P.T.MD,P.T.MD+P.T.recovery]); xlim([-P.T.baseline,P.T.MD+P.T.recovery]); grid on;
xlabel('time (day) '); ylabel('y'); ylim([.5,1.5]);

%%% Functions %%%
function CopyStruct(In)
% copy fileds of a structure array to local variables
flds=fieldnames(In);
for i=1:length(flds)
    assignin('caller', flds{i}, In.(flds{i}));
end

function Out=FinalValues(In)
% extract the final part of structure
flds=fieldnames(In);
for i=1:length(flds)
    Out.(flds{i})=In.(flds{i})(:,end);
end

function Out=AppendStruct(In1,In2)
% append In2 to In1
flds=fieldnames(In1);
for i=1:length(flds)
    Out.(flds{i})=[In1.(flds{i}), In2.(flds{i})];
end

function Out=OutSub(In1,In2)
% outer product of the subtraction of two vectors 
if size(In1,1)<size(In1,2); In1=In1'; end;
if size(In2,1)<size(In2,2); In2=In2'; end;
Out=In1*ones(size(In2'))-ones(size(In1))*In2';

function [Mu,Q,S]=InputStat(N,Input,MDflag)
% generate input mean vector and covariance matrix
CopyStruct(Input);
F=@(x) 1./(1+exp(3*(x-1)));
muC=mu*(MD*(MDflag==1)+1*(MDflag~=1));
muI=mu*(MD*(MDflag==2)+1*(MDflag~=2));
axC=transpose(linspace(0,1,N.contra));
axI=transpose(linspace(0,1,N.ipsi));
Mu=[muC*ones(size(axC)); muI*(ones(size(axI)))];
QCC=muC^2*exp(-.5*OutSub(axC,axC).^2/L^2)/sqrt(2*pi*L^2);
QII=muI^2*exp(-.5*OutSub(axI,axI).^2/L^2)/sqrt(2*pi*L^2);
QCI=BE*muC*muI*exp(-.5*OutSub(axC,axI).^2/L^2)/sqrt(2*pi*L^2);
Q=[QCC,QCI; QCI',QII];
S=F([(axC-.5).^2; (axI-.5).^2]./L^2);
S=S./sum(S);
Nz=randn(N.total,N.total); 
Q=Q+(Nz+Nz'); 

function Hist=PlastSim(V,S,Flags,P,stage)
% simulation the two-factor rule of plasticity
TL=@(x) x.*(x>0);
N=P.N.total;
X0=[V.rho; V.Hh]; %initialization
[Mu,Q]=InputStat(P.N,P.Input,Flags.MD); %acquire input statistics
[T,X]=ode45(@(t,X) PlastRule(t,X,Flags,P,Mu,Q,S),V.time(end)+[0,P.T.(stage)],X0);
X=transpose(X);

Hist.time=transpose(T(2:end));
Hist.rho=X(1:N,2:end); 
Hist.Hh=X(end,2:end);
Hist.H=P.Hfun(X(end,2:end));
if Flags.TNF==0; Hist.H=ones(size(X(end,2:end))); end;
Hist.w=Hist.rho.*(ones(N,1)*Hist.H).*(S*ones(size(Hist.time)));
Hist.y=Mu'*Hist.w;
Hist.resp=[sum(Hist.w(1:P.N.contra,:)); sum(Hist.w(P.N.contra+1:P.N.total,:))];
Hist.ODI=(Hist.resp(1,:)-Hist.resp(2,:))./(Hist.resp(1,:)+Hist.resp(2,:));

function dX=PlastRule(t,X,Flags,P,Mu,Q,S)
% dynamical equations used in PlastSim
TL = @(x) x.*(x>0); %threshold-linear function
N=P.N.total;
rho=X(1:N);
Hh=X(end);
H=P.Hfun(Hh);
if Flags.TNF==0; H=1; end;
w=H*rho.*S;
y=Mu'*w;
Corr=Q*w-P.th.corr;

HON=(Flags.NMDA)||(t<4);
dX=zeros(N+1,1);
dX(1:N)=HON*(Flags.Trk*TL(P.rhoM-rho).*TL(Corr)-TL(rho-P.rhom./sqrt(H)).*TL(-Corr))/P.tau.rho;
dX(end)=(-Hh+P.Hhfun(H*P.y0/y))/P.tau.Hh;
