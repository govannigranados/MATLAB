function [y,Time,yPrevSub,stats]= TimeAdvance(y,t,options,F,...
                                  method,SubStep,h,Time,Cycles,yPrevSub)
    stats=0;                                
    switch(method)
    %------------------------------------------
    %Begin time integration options
    %------------------------------------------
    
    
    %------------------------------------------
    %Built in options
    %------------------------------------------
    
    %------------------------------------------
    %ODE45          Opt=45
    %------------------------------------------
        case('45')
            for j=1:Cycles
                [tt,U]=ode45(F, [y(end),y(end)+SubStep],y,options);
                y=U(end,:);
                Time=[Time;tt(2:end)];
            end
    %------------------------------------------
    %ODE15s Opt=15s, Status: working
    %------------------------------------------
        case('15s')
            for j=1:Cycles %Cycles for memory reasons
                [tt,U]=ode15s(F, [y(end),y(end)+SubStep],y,options);
                y=U(end,:);
                Time=[Time;tt(2:end)];
            end
            
     
    %------------------------------------------
    %Implicit options
    %------------------------------------------
           
    %------------------------------------------
    %Adams Moulton 1 Opt=AM1, Status:Broken
    %------------------------------------------
        case('AM1')
            for j=1:Cycles
                G=@(x)(y-x+SubStep*F(t,x));%Time independent, shenangins
                guess=y;
                options = optimoptions('fsolve','Display','none');
                [y,fval]=fsolve(G,guess, options);
            end
    %------------------------------------------
    %Adams Moulton 2 Opt=AM2, Status:Uncertain
    %------------------------------------------
        case('AM2')
            G=@(x)(y-x+h/2*(F(t,x)+F(t,y)));%Time independent, shenangins
            guess=y;
            options=optimoptions('fsolve', 'FunctionTolerance', 1e-2...
            ,'StepTolerance', 1e-2,'OptimalityTolerance', 1e-2);
            [y,fval]=fsolve(G,guess, options);
    
    %------------------------------------------
    %BDF2 2 Opt=BDF2, Status: Tentatively working,
    %------------------------------------------
        case('BDF2')
            
            options=optimoptions('fsolve','Display','none','FunctionTolerance',1e-8,...
                        'StepTolerance', 1e-8,'OptimalityTolerance', 1e-8);
            %y_(n+2)-4/3 y_(n+1)+1/3 y_n = 2/3 f(t_(n+2),y_(n+2))
            for j=1:Cycles
                if(yPrevSub==y)    
                    G=@(x) x-y-SubStep*F(t,x);
                else
                    G=@(x) x-4/3*y+1/3*yPrevSub-2/3* SubStep*F(t,x);
                end 
                %}
                [yNew,~]=fsolve(G,y,options);
                yPrevSub=y;
                y=yNew;
            end
            %Needs a solver and previous step
    %------------------------------------------
    %NewtonKrylovBisetti Opt=NKB, Status: In progress, reduced to first order
    %------------------------------------------
        case('NKB')
            if stats ==0
            clear stats
            end
            %y_(n+2)-4/3 y_(n+1)+1/3 y_n = 2/3*f(t_(n+2),y_(n+2))
             warning('off','all')
             opts.newton.tol = 1e-9;%try 1e-9 on advancing 1e3
             opts.newton.d   = 1; 
             opts.newton.maxiters = length(y); 
             opts.newton.inexact =1;
             opts.gmres.restart=length(y);
             opts.gmres.maxiters=length(y)^2;
             opts.gmres.tol=1e-8;
             opts.gmres.tolbnds = [1e-9,1e-7];%[1e-12,1e-9]
            for j=1:Cycles
                if(yPrevSub==y)    
                    G=@(t,x) x-y'-SubStep*F(t,x);
                else
                    G=@(t,x) x-4/3*y'+1/3*yPrevSub'-2/3* SubStep*F(t,x);
                end   
                [yNew,OutStats]=NewtonKrylovPGMRES(y(end),y',G,[],opts);
                yPrevSub=y;
                y=yNew';
                Time=[Time,y(end)];
                stats.newtonIter(j)=OutStats.newton.nriters;
                stats.gmresIter(j,:)=OutStats.gmres.iter;
            end
            %Needs a solver and previous step
            
            
    %------------------------------------------
    %NewtonKrylov Opt=NK, Status: In progress, reduced to first order
    %------------------------------------------
        case('NK')
            %y_(n+2)-4/3 y_(n+1)+1/3 y_n = 2/3*f(t_(n+2),y_(n+2))
            for j=1:Cycles
                if(yPrevSub==y)    
                    G=@(t,x) x-y-SubStep*F(t,x);
                else
                    G=@(t,x) x-4/3*y+1/3*yPrevSub-2/3* SubStep*F(t,x);
                end   
                [~,JacG]=ApproxJac(y,t,G,0);
                [yNew,~]=gmres(JacG,JacG*y'-G(t,y)',length(JacG),1e-4);
                yPrevSub=y;
                y=yNew';
                Time=[Time,y(end)];
            end
            %Needs a solver and previous step
    
            
    %------------------------------------------
    %Explicit methods
    %------------------------------------------
    
    %------------------------------------------
    %Backward Euler Opt=Euler,
    %------------------------------------------
        case('Euler')
        y=y+h*F(t,y);
    %------------------------------------------
    %RK2 Opt=RK2  
    %------------------------------------------
        case('RK2')
        Time=y(end);
        for j=1:Cycles    
            y=y+SubStep*F(t+.5*SubStep,y+.5*SubStep*F(t,y));
            Time=[Time,y(end)];
        end
    %------------------------------------------
    %RK4  Opt=RK4           Status: Working
    %------------------------------------------
        case('RK4')
            Time=y(end);
            for j=1:Cycles
                k1=SubStep*F(t,y);
                k2=SubStep*F(t+h/2,y+k1./2);
                k3=SubStep*F(t+h/2,y+k2./2);
                k4=SubStep*F(t+h,y+k3);
                y=y+(1/6)*(k1+2*k2+2*k3+k4);
                Time=[Time,y(end)];
            end
            
    %------------------------------------------
    %End time integration options
    %------------------------------------------        
    end

end

