function [F,X] =Load_ShiftingV3(X, Pon, Poff, Tonmax)
 F=0;
 Nc=length(X(:,1));
 Ts=length(X(1,:));
 P_total=Pon*X+Poff'*(1-X);
  peak=find(P_total>=max(P_total)); 
  p=1;
  n=1;
 while(p<=length(peak))
   while (n<=Nc)
        if(X(n,peak(p))==1)
            Gate=zeros(Ts,1);
            Ton_Gate=max(1,peak(p)-Tonmax+1):min(peak(p)+Tonmax-1,Ts);
            Gate(Ton_Gate)=1;
            X1=Gate.*X(n,:)';
            shift=fix(sum(X1)/2)+1; % Need to think
            X2=zeros(size(X1));
            X2(1:end-shift)=X1(shift+1:end); 
            X_new=X;
            X_new(n,:)=X2+((1-Gate).*X(n,:)');
           P_total_new=Pon*X_new+Poff'*(1-X_new);
           if( max(P_total_new)<max(P_total)&&sum(Gate([1:shift]))==0)
               X=X_new;
               F=1;
               n=Nc+1;
               p=length(peak)+1;
               break
           end
           
           X2(shift+1:end)=X1(1:end-shift); 
           X_new=X;
           X_new(n,:)=X2+((1-Gate).*X(n,:)');
           P_total_new=Pon*X_new+Poff'*(1-X_new);
           if( max(P_total_new)<max(P_total))
               X=X_new;
               F=1;
               n=Nc+1;
               p=length(peak)+1;
               break
           end
          
        end
         n=n+1;
   end
     p=p+1;
 end

 
 end