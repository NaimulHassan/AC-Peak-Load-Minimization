clc
close all
clear all
for iter=1
tic
Nc=100;% nubmber of AC
Tout=32*ones(Nc,1);
Tset=zeros(Nc,1);
A=zeros(Nc,1);
Tset_Set=16:28;
for n=1:Nc
    Tset(n)=Tset_Set(ceil(rand*length(Tset_Set)));
end
%load(['Hybrid\' num2str(Nc) '\archive\input\Tset_' num2str(iter) '.mat']);

load(['Tset' '.mat']);

delT_low=2;
delT_high=6;
Area_Set=[30 40 50 55]; %[37.16 51.0967 111.49 167.7849];
C_Set= [1 1.5 2 3];


for n=1:Nc
    a=ceil(rand*length(Area_Set));
    A(n)=Area_Set(a);
    C(n)=C_Set(a);
end

%load(['Hybrid\' num2str(Nc) '\archive\input\C_' num2str(iter) '.mat']);
load(['C' '.mat']);
for(i=1:length(C))
    if C(i)==3
        A(i)=55;
    elseif C(i)==2
        A(i)=50;
    elseif C(i)==1.5
        A(i)=40;
    else
        A(i)=30
    end
end


%save(['Data Extraction\' num2str(Nc) '\input\C_' num2str(iter) '.mat'], 'C');
COP=2.9;

%% Constant values
Req=0.35;
Cp=1.005*10^3;
Density_Air=1.225;
h= 3.2;
v= A*h;
Mair= Density_Air*v;eq= 0.35;
Ph_on= -3517*C;
Ph_off= 373*COP*ones(Nc,1);
Pon= ceil(-Ph_on/COP);
Poff= Ph_off/COP;

TON_high= (-Cp.*Mair.*Req.*delT_high./(Ph_on'.*Req+Tout-Tset))./60;
TOFF_high= abs(Cp.*Mair.*Req.*delT_high./(Ph_off.*Req+Tout-Tset))./60;

TON_low= (-Cp.*Mair.*Req.*delT_low./(Ph_on'.*Req+Tout-Tset))./60;
TOFF_low= abs(Cp.*Mair.*Req.*delT_low./(Ph_off.*Req+Tout-Tset))./60;

%%%%% Here save the values for the instance 
% save intance1 Nc, Tset, Tout, A, C;

%%%%%%%%%%%%%  start the heuristic scheduling aldorithm %%%%%%%%%%%%%

Ts=90;% slot numbers
Tonmax=ceil(TON_high);
Tonmin=ceil(TON_low); 
Toffmax=fix(TOFF_high);
Toffmin=fix(TOFF_low);

C_Set=sort(C_Set,'descend');
[~,s]=sort(C,'descend');
C=C(s);
TON_high=TON_high(s);
TOFF_high=TOFF_high(s);
TON_low=TON_low(s);
TOFF_low=TOFF_low(s);
Tonmax=Tonmax(s);
Tonmin=Tonmin(s);
Toffmax=Toffmax(s);
Toffmin=Toffmin(s);

K=length(C_Set);% types of AC
N=zeros(K,Nc);% No of Acs in type k

for n=1:Nc
        for k=1:K
            if(C(n)==C_Set(k))
                N(k,n)=1;
                            
              end
        
    end
end
N_r=N;% left for scheduling
X=zeros(Nc,Ts);

%%%%%%%%%%%%%% Scheduling of same type of ACs one after another
Nss=fix((TON_high+TOFF_high)./TON_high);
for k=1:K
    Ns(k)=max(Nss.*N(k,:)');% no of 
    Ton_s(k)=max(Tonmax.*N(k,:)');
end

Toff_s=(Ns-1).*Ton_s;
Tc_s=Ton_s+Toff_s;

for k1=1:K
      M=0;
      if(sum(N(k1,:))>=Ns(k1))
      %check=1
      M=fix(sum(N(k1,:))/Ns(k1));
       u= find(N(k1,:)==1);
      for m=1:M           
         for ns=1:Ns(k1)
           %ns
           %k1
           t_on=Ton_s(k1)*(ns-1)+[1:Ton_s(k1)];
           X(u(ns),t_on)=1;
           N_r(k1,u(ns))=0;
           count1=1;
            while(max(t_on)<Ts)
                %check=2
                 t_on=Ton_s(k1)*(ns-1)+[1:Ton_s(k1)]+count1*Tc_s(k1);
                  count1=count1+1;
                 if(max(t_on)<=Ts)
                       X(u(ns),t_on)=1;
                 end
                 if(max(t_on)>Ts)
                     for t1=1:length(t_on)
                       if (t_on(t1)<=Ts) 
                           X(u(ns),t_on(t1))=1; 
                       end
                     end
                 end
            end
             
         end
      end
      
      end
          
 end

%%%%%%%%%%%%%% Scheduling the AC with higher capacity one after another 
  k_r=1;
  while(sum(N_r(k_r,:))==0)
      k_r=k_r+1;
  end
check_Nr=N_r;
check_k=k_r ;
if(sum(N_r(k_r,:))>0)
  u1= find(N_r(k_r,:)==1) ;    
         for ns1=1:sum(N_r(k_r,:))
           t_on=Ton_s(1)*(ns1-1)+[1:Ton_s(1)];
           X(u1(ns1),t_on)=1;
           N_r(k_r,u1(ns1))=0;
           count2=1;
            while(max(t_on)<Ts)                
                 t_on=Ton_s(1)*(ns1-1)+[1:Ton_s(1)]+count2*Tc_s(1);
                  count2=count2+1;
                 if(max(t_on)<=Ts)
                       X(u1(ns1),t_on)=1;
                 end
                 if(max(t_on)>Ts)
                     for t2=1:length(t_on)
                       if (t_on(t2)<=Ts) 
                           X(u1(ns1),t_on(t2))=1; 
                       end
                     end
                 end
            end
         end
end
    
P_total=(sum(N-N_r).*Pon)*X+(sum(N-N_r).*Poff')*(1-X);
    
   %%%%%%%%%%%%%% Scheduling of combination of different types of ACs one after another
   del_ON=Tonmax-Tonmin;
   
 for k2=k_r+1:K
     if(sum(N_r(k2,:))>0)
         u3=find(N_r(k2,:)==1);
         for n3=1:length(u3)
             t3=1;
             t_delay=0;
             last_Ton_Time=1;
             while(t3<Ts)
                 true=zeros(del_ON(u3(n3))+1,1);
                 for dt=1:del_ON(u3(n3))+1
%                      n3
%                      u3
%                      k2
                     %X
                     a=t3:min(Ts,t3+Tonmin(u3(n3))+dt-2);
                     b=min(Ts,t3+Tonmin(u3(n3))+dt-1):min(Ts,t3+Tonmin(u3(n3))+dt-1+(Ns(k2)-1)*(Tonmin(u3(n3))+dt-1));
                    if(max(P_total(a))<= min(P_total(b)))
                     true(dt)=1;
                     %true
                    end
                 end
                 
                 
%                  t3
%                  sum(true)
                 if(sum(true)>=1)
                   flag=2;
                   %P_total
                   t_on=t3:min(Ts,t3+Tonmin(u3(n3))+sum(true)-2);
                   X(u3(n3),t_on)=1;
                   %P_total=(sum(N-N_r).*Pon)*X+(sum(N-N_r).*Poff')*(1-X)
                   t3=t3+length(t_on)*Ns(k2);
                   k2;
                   t_delay=length(t_on)*(Ns(k2)-1);
                   last_Ton_Time=length(t_on);
                 end
                  
                   if(sum(true)<1)
                         t_delay;
                         t3;
                         Ns;
                         k2;
                         c=t3:min(Ts,t3+fix((t_delay)/(Ns(k2)-1))-1);
                         d=min(Ts,t3+fix((t_delay)/(Ns(k2)-1))):min(Ts,t3+fix((t_delay)/(Ns(k2)-1))+t_delay);
                      if(t_delay>(Ns(k2)-1)*Tonmin(u3(n3))&& max(P_total(c))<max(P_total(d)))
                         flag=1;
                         t_on=t3:t3+min(Tonmax(u3(n3)),fix((t_delay)/(Ns(k2)-1)))-1;                  
                          X(u3(n3),t_on)=1;
                          t3=t3+length(t_on)*Ns(k2)-1;
                          t_delay=length(t_on)*(Ns(k2)-1)-1;
                          last_Ton_Time=length(t_on);
                       end     
                       
                        if(t_delay>(Ns(k2)-1)*Tonmax(u3(n3)))
                         flag=0;
                         t_on=t3:t3+Tonmax(u3(n3))-1          ;             
                          X(u3(n3),t_on)=1;
                          t3=t3+length(t_on)*Ns(k2)-1;
                          t_delay=length(t_on)*(Ns(k2)-1)-1;
                          last_Ton_Time=length(t_on);
                        end 
                       
                   t3=t3+1; 
                   t_delay=t_delay+1;
                 end
             end
              N_r(k2,u3(n3))=0;
              P_total=(sum(N-N_r).*Pon)*X+(sum(N-N_r).*Poff')*(1-X);
         end
     end
 end

 max(P_total);
 %%%%%%%%%%%%%%%%%% shifting the loads of peak point
  F=1;
 while(F==1)
 %[F,X] =Load_Shifting(X, Pon, Poff, Tonmax);  % change the function here
 [F,X] =Load_ShiftingV3(X, Pon, Poff, Tonmax);
 end
execution_time_heuristic(iter)= toc;
P_total=Pon*X+Poff'*(1-X);
P_total'
Total_P_heuristic(iter)= sum(P_total/1000);
Max_P_heuristic(iter)= max(P_total/1000)
Var_P_heuristic(iter)= var(P_total/1000);

end

mean_execution_time= mean(execution_time_heuristic);
avg_total_power= mean(Total_P_heuristic);
avg_max_power= mean(Max_P_heuristic)
avg_var= mean(Var_P_heuristic)

Max_P_heuristic';
%S_M_P= var(Max_P_heuristic);
%S_M_P/avg_max_power


%%save
% save(['E:\Academic_Life\4-2\Thesis\Optimization Problem\Sir files\results' '\Total_P_heuristic_' num2str(Nc) '.mat'], 'Total_P_heuristic')
% save(['E:\Academic_Life\4-2\Thesis\Optimization Problem\Sir files\results' '\Max_P_heuristic_' num2str(Nc) '.mat'], 'Max_P_heuristic')
% save(['E:\Academic_Life\4-2\Thesis\Optimization Problem\Sir files\results' '\Var_P_heuristic_' num2str(Nc) '.mat'], 'Var_P_heuristic')
% save(['E:\Academic_Life\4-2\Thesis\Optimization Problem\Sir files\results' '\execution_time_heuristic_' num2str(Nc) '.mat'], 'execution_time_heuristic')