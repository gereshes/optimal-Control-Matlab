close all
clear all
clc

%% Variable Setup
 
 
 
%In km
m0=3366;
rEarth= 6378.14;%km
mu=398600.44;%(km**3 s**-2)
 
 
a=rEarth+400 %For the ISS
 
%Setting up characteristic quantities
lStar=a
tStar=sqrt((lStar^3)/mu);
mStar=m0
 
 
c1=1.33/(1000); %Thrust in Newtons (using 7RIT2X engine)
c1=c1*tStar*tStar/(lStar*mStar); %Non-dimensionalizing the problem
isp=3872; %seconds(using RIT2X engine)
g0=9.81/1000; %Must be in km s**-2 to be consistent with the rest of the constants 
c2=(isp*g0);
c2=c2*tStar/lStar;

n=sqrt(mu/(a^3))*tStar


args=[c1,c2,n];

entries=65%500
t0=0
%tof=1000/tStar
tof=1000/tStar
d=1
xTarget=[1,0,0,1]

%xChaser=[1.1,-.01,0,1.01*sqrt(1/1.1)]
xChaser=[.9,0,0,sqrt(1/.9)]
%xChaser=[1.1,-.01,-.0001,1.1*sqrt(1/1.1)]
x0=[xTarget-xChaser,1]
numStates=5;
numControls=2;
%x0=[.1,0,0,.1,1]

%% Forming the Initial Guess

U=ones(2,entries);
U(1,:)=U(1,:)*-pi;
U(2,:)=U(2,:)*.4;

tAct=linspace(0,tof,entries);

%sol=ode45(@(t,states) dynamicsTemplateInitGuess(t,states,args),[t0,t0+tof],x0);

xf=[0/lStar,0/lStar,0*tStar/lStar,0*tStar/lStar,.1]
setConds=[x0,t0;xf,t0+tof];

%X=[deval(sol,tAct);tAct];
%X(1:5,end)=xf';

X=[linspace(x0(1),xf(1),entries);...
   linspace(x0(2),xf(2),entries);...
   linspace(x0(3),xf(3),entries);...
   linspace(x0(4),xf(4),entries);...
   linspace(x0(5),xf(5),entries);...
   tAct];
%X=[1.52305651741452,1.52310509022741,1.52295785829696,1.52212128092207,1.52015765837415,1.51666184927679,1.51125711464935,1.50356520108333,1.49337792755994,1.48046463816510,1.46460278321186,1.44577036479365,1.42394907141009,1.39915240701545,1.37146925821388,1.34126512701403,1.30891267583629,1.27485323112760,1.23964541897514,1.20387423382040,1.16821448530000,1.13335503967040,1.10025906656605,1.07026892976355,1.04512551059061,1.02578526652402,1.01237991271947,1.00439667507896,1.00079724442264,0.999999999999973;-4.48836528787399e-12,0.126124247135746,0.251227045862207,0.375390003820761,0.498754490214382,0.621549835944247,0.744078268895695,0.866711872418154,0.989889933262621,1.11408596920300,1.23986429806745,1.36783437873105,1.49867797372117,1.63313228625189,1.77200595362509,1.91618945178296,2.06656111410582,2.22407821801673,2.38966046571926,2.56416632659858,2.74869945126129,2.94430296442545,3.15110149937880,3.36951042585093,3.59811763890243,3.83448429625227,4.07592550098173,4.31951253178059,4.56243978519018,4.80253374041831;-4.10560055938904e-12,0.000125909961585124,-0.00173556619170321,-0.00558800879584833,-0.0112008438306392,-0.0184183936201109,-0.0273070118894244,-0.0374261340692132,-0.0483398576903365,-0.0603682985015653,-0.0728804117309810,-0.0853955074964949,-0.0979423167642875,-0.110411905286094,-0.121880500324341,-0.131678928667053,-0.139869950538973,-0.145906684854686,-0.149560507382840,-0.150496021069847,-0.148667948084767,-0.143504502104040,-0.133821586329831,-0.116619142461490,-0.0938867973803568,-0.0686091633657524,-0.0443980362871125,-0.0233889265227163,-0.00799714431683030,-7.42694798164085e-15;0.810272181148086,0.803801239133021,0.797337914233086,0.791277396719283,0.785792361331415,0.781256775915317,0.777870774222825,0.775926377213847,0.775503180660884,0.776802949557059,0.780043356344445,0.785323792832707,0.792815932005139,0.802515092407948,0.814645427959141,0.829071037643674,0.845698065599929,0.864479429742497,0.884942996396127,0.906419413392498,0.932924230155431,0.957887118754226,0.982100167087131,1.00753128553441,1.02252552837196,1.03171293305918,1.03360265743820,1.02832996646797,1.01673469886311,1.00059168177094;0.999999999958119,0.994906839421100,0.989795757908893,0.984764916118786,0.979757240610204,0.974897852471694,0.970172400964529,0.965660976037318,0.961343288446346,0.957208758100501,0.953321647968988,0.949699386702377,0.946346553805600,0.943145648608485,0.940153742968048,0.937193381133902,0.934195358669008,0.931113740408267,0.927797775048685,0.923986625694904,0.922877849431475,0.925544226216827,0.928839169238942,0.927043129366068,0.919924019035006,0.912018222100821,0.903181653675787,0.893304679133800,0.882348453693619,0.870436701192429]


tSpan=linspace(0,tof,100);
xPosIG=interp1(tAct,X(1,:),tSpan,'spline');
heightIG=interp1(tAct,X(2,:),tSpan,'spline');





%% Enforcing forward Time and intial and final states
statesPerEntry=entries*(numStates+numControls+1);
A=zeros(entries,statesPerEntry);
for c=1:entries-1
    for d=1:statesPerEntry
        if(d==(c*(numStates+1)))
            A(c,d)=1;
        elseif(d==((c+1)*(numStates+1)))
            A(c,d)=-1;
        end
    end
end

A=sparse(A);
B=-ones(entries,1)*eps;
Aeq=zeros(5+1+4,statesPerEntry);
for c=1:6 %Enfocing intial state and time
    Aeq(c,c)=1;
end
for c=1:4 %enforcing final state
    Aeq(6+c,((entries-1)*(numStates+1))+c)=1;
end
Beq=[x0';t0;xf(1:4)'];
linear={A,B,sparse(Aeq),Beq};

%% Forming the Bounds
uBDS=[Inf,Inf,Inf,Inf,1.01,Inf]
lBDS=[-Inf,-Inf,-Inf,-Inf,.01,-.1]
uBDC=[4*pi,1]
lBDC=[-4*pi,1e-6]

bds={repmat(lBDS',1,entries),repmat(uBDS',1,entries),repmat(lBDC',1,entries),repmat(uBDC',1,entries)}
fun=@(x,u,args)dynamicsTemplate(x,u,args)

%% Solving the Optimal Control Problem
[X,U,tAct, Fval, ExitFlag, Output]=fminconWrapper(X,U,fun,args,bds,linear);
%assert(ExitFlag==1, 'The solver did not find a feasable solution')


%% Plotting the solution
control=[U;tAct];
tSpan=linspace(0,tAct(end),100);
%tAct=linspace(0,tof,entries)
xRelPos=interp1(tAct,X(1,:),tSpan,'spline');
yRelPos=interp1(tAct,X(2,:),tSpan,'spline');
xDotRel=interp1(tAct,X(3,:),tSpan,'spline');
yDotRel=interp1(tAct,X(4,:),tSpan,'spline');
mRel=interp1(tAct,X(5,:),tSpan,'spline');
thrust=interp1(tAct,U(2,:),tSpan,'linear');
angle=interp1(tAct,U(1,:),tSpan,'linear');

figure()
hold on
plot(xRelPos,yRelPos,'k','DisplayName','Trajecotry')

figure()
plot(tSpan,thrust)
hold on
plot(tSpan,angle)





%% Investigating Error
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
solComp=ode45(@(t,states) dynamicsTemplateForComp(t,states,args,control),tSpan,x0,opts);

compStates=deval(solComp,tSpan);

errorStates=compStates-[xRelPos;yRelPos;xDotRel;yDotRel;mRel];
vecNorm=vecnorm(errorStates)
figure()
plot(tSpan,vecNorm)
set(gca, 'YScale', 'log')