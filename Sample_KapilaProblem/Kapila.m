function [y,TIME]=Kapila(a,b,e,l0,Tend,Delt,method)
%e is epsilon
L=l0*exp(1/e);
z=e*b;
z=5e-3;
y0=[1,a-z,z]; % theta(0),y(0),z(0)
y=y0;
if method=="15s"
    F=@(t,y)[rhs(t,y,e,L)]';
elseif method=="45"
    F=@(t,y)[rhs(t,y,e,L)]';
else
    F=@(t,y)[rhs(t,y,e,L)];
end
NumSteps=round(Tend/Delt);
options=odeset('RelTol', 1e-18,'AbsTol',1e-18);
Cycles=1;
h=Delt;
SubStep=Delt;
Time=0;
yPrevSub=y;
start=tic;
%Start the substepping
for i=1:NumSteps
    t=y(end);
    [yUp,Time,yPrevSub,stats]= TimeAdvance(y(end,:),t,options,F,...
                                  method,SubStep,h,Time,Cycles,yPrevSub);
%Clean data
    for j=1:3
        if(yUp(j)<0)
            yUp(j)=0;
        end
    end
    y(i+1,:)=yUp;
    fprintf(['Completed time step:', num2str( (i)*h), '\n'])
end
TIME=toc(start);
y=y';
time=0:h:10;
figure(1)
plot(time,y,"LineWidth",1.5,"MarkerSize",15)
xlabel("Time")
ylabel("Non-Dimensional quantity")
title('states over time')
legend({'Temp','y','z'});

for i=1:length(y)
   J=jac(y(end,i),y(1:3,i),e,L);
   eig_(i,:)=sort(eig(J));
end

figure(2)
semilogy(time,abs(real(eig_)),"LineWidth",1.5,"MarkerSize",15);
title('real part of eigs over time')
xlabel('time')
ylabel('real part')
legend({'E1','E2','E3'});


function yd=rhs(t,y,e,L)

yd=zeros(size(y)-1);

omega=L*y(2).*y(3).*exp(-1./(e*y(1)));
yd(1)=y(3);
yd(2)=-omega;
yd(3)=-yd(2)-yd(1);


% y = [theta,y,z]
function J=jac(~,y,e,L)

n=length(y);
J=zeros(n,n);

%omega=l*y(2).*y(3).*exp(-1./(e*y(1)));

J(1,3)=1;
J(2,1)=-(L*y(2)*y(3)*exp(-1/(e*y(1))))/(e*y(1)^2);
J(2,2)=-L*y(3)*exp(-1/(e*y(1)));
J(2,3)=-L*y(2)*exp(-1/(e*y(1)));
J(3,1)=-J(2,1);
J(3,2)=-J(2,2);
J(3,3)=-J(2,3)-1;



