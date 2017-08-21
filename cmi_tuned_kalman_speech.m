function [cleanSpeech] = cmi_tuned_kalman_speech(x,fs)

%Algorithm as described in the paper, "Application of the tuned kalman 
%filter in speech enhancement." IEEE CMI 2016 - O.Das, B Goswami and R
%Ghosh

%Inputs:
%x - noisy speech sample
%fs - sampling rate
%Output:
%cleanSpeech - enhanced speech

p=15;%order of lpc
y=x';

%dividing into overlapping 80ms frames
start=1;
l=0.08*fs;
overlap=0.01*fs;
totseg=ceil(length(y)/(l-overlap));
segment=zeros(totseg,l);

for i=1:totseg-1
    segment(i,1:l)=y(1,start:start+l-1);
    start=(l-overlap)*i+1;
end
segment(totseg,1:length(y)-start+1)=y(start:length(y));


H=[zeros(1,p-1),1];
G=H';
cleanspeech=zeros(totseg,l);
cleanSpeech=zeros(1,length(y));


[R, silent_inds] = measurementNoiseNew(segment, fs);       
J1=zeros(1,10);
J2=zeros(1,10);
nq=zeros(1,10);
u=1;
X=y(1:p)';
P=zeros(l,p,p);
P(1,:,:)=R*eye(p);
Q=0;
t1=zeros(p,p);

K_Q2=zeros(p,l,totseg);
K_Q3=zeros(p,l,totseg);
Kavg_Q2=zeros(1,totseg);
Kavg_Q3=zeros(1,totseg);
J1_Q3=zeros(1,totseg);
J1_Q2=zeros(1,totseg);
Xhat_Q2=zeros(p,l,totseg);
Xhat_Q3=zeros(p,l,totseg);


for i=1:totseg
    
    %first iteration of Kalman filter
    [A,Q1]=lpc(segment(i,:),p);
    temp=eye(p);
    PHI=[temp(2:p,:);-fliplr(A(2:end))];
    
    %tuning the filter by calculating optimum value of process noise
    %variance
   q=1;
   if i~=1

     for n=-5:4
       Q0=(10^n)*Q1;
       t1(:,:)=P(1,:,:);
       Ak=H*(PHI*t1*PHI')*H';
       Bk=H*Q0*H';
       J1(q)=R/(Ak+Bk+R);
       J2(q)=Bk/(Ak+Bk);
       nq(q)=log10(Bk);
       q=q+1;
     end
     
     
     [nq_nom,J1_Q2(i)]=intersections(nq,J1,nq,J2);
     
     %sensitivity metrics gets higher preference and we select Q<Qcomp
     if numel(nq_nom)~=0
            Q=10^(nq_nom-0.7);
     else
       Q=Q1;
     end
        
   else
     Q=Q1;
   end
   u=u+1;
   
   for j=1:length(segment(i,:))
        X_=PHI*X;
        t1(:,:)=P(j,:,:);
        P_=(PHI*t1*PHI')+(G*Q*G');
        K_Q2(:,j,i)=(P_*H')*(inv(H*P_*H'+R));
        t1=(temp-K_Q2(:,j,i)*H)*P_;
        P(j+1,:,:)=t1(:,:);
        e=segment(i,j)-(H*X_);
        X=X_+K_Q2(:,j,i)*e;
        cleanspeech(i,j)=X(end);
        Xhat_Q2(:,j,i)=X_;
        
   end
    Kavg_Q2(i)=mean(K_Q2(p,:,i));
    P(1,:,:)=P(j-1,:,:);
       
    
end
u=1;

 for i=1:totseg
    
    %first iteration of Kalman filter
    [A,Q1]=lpc(segment(i,:),p);
    temp=eye(p);
    PHI=[temp(2:p,:);-fliplr(A(2:end))];
    
    %tuning the filter by calculating optimum value of process noise
    %variance
   q=1;
   if i~=1

     for n=-5:4
       Q0=(10^n)*Q1;
       t1(:,:)=P(1,:,:);
       Ak=H*(PHI*t1*PHI')*H';
       Bk=H*Q0*H';
       J1(q)=R/(Ak+Bk+R);
       J2(q)=Bk/(Ak+Bk);
       nq(q)=log10(Bk);
       q=q+1;
     end
     
     
     [nq_nom,J1_Q3(i)]=intersections(nq,J1,nq,J2);
     
     %sensitivity metrics gets higher preference and we select Q<Qcomp
     if numel(nq_nom)~=0
            Q=10^(nq_nom);
     else
       Q=Q1;
     end
        
   else
     Q=Q1;
   end
   u=u+1;
   
   for j=1:length(segment(i,:))
        X_=PHI*X;
        t1(:,:)=P(j,:,:);
        P_=(PHI*t1*PHI')+(G*Q*G');
        K_Q3(:,j,i)=(P_*H')*(inv(H*P_*H'+R));
        t1=(temp-K_Q3(:,j,i)*H)*P_;
        P(j+1,:,:)=t1(:,:);
        e=segment(i,j)-(H*X_);
        X=X_+K_Q3(:,j,i)*e;
        cleanspeech(i,j)=X(end);
        Xhat_Q3(:,j,i)=X_;
        
   end
    Kavg_Q3(i)=mean(K_Q3(p,:,i));
    P(1,:,:)=P(j-1,:,:);
        
 end

K=zeros(1,totseg);
I=ones(p,1);

for i=1:totseg
         
      if(silent_inds(i) == 1)
       %silent region
          for j=1:length(segment(i,:))
              X=K_Q2(:,j,i)*segment(i,j)+(I-K_Q2(:,j,i))*(H*Xhat_Q2(:,j,i));
              cleanspeech(i,j)=X(end);
          end
          K(i)=Kavg_Q2(i);
      else
      %voiced region
          for j=1:length(segment(i,:))
              X=K_Q3(:,j,i)*segment(i,j)+(I-K_Q3(:,j,i))*(H*Xhat_Q3(:,j,i));
              cleanspeech(i,j)=X(end);
          end
           K(i)=Kavg_Q3(i);
      end
     
     %second iteration of kalman filter with lpcs calculated from cleaned speech      
     [A,Q]=lpc(cleanspeech(i,:),p);
     temp=eye(p);
     PHI=[temp(2:p,:);-fliplr(A(2:end))];
     if i==1
        X=cleanspeech(i,1:p)';
        P0=temp*R;
     end
    
    for j=1:length(segment(i,:))
        
        X_=PHI*X;
        P_=(PHI*P0*PHI')+(G*Q*G');
        K0=(P_*H')*(inv(H*P_*H'+R));
        P0=(eye(p)-K0*H)*P_;
        e=segment(i,j)-(H*X_);
        X=X_+K0*e;
        cleanspeech(i,j)=X(end);
    end
       
end



%overlap add
cleanSpeech(1:l)=cleanspeech(1,1:l);
start=l+1;
for i=2:totseg-1
    cleanSpeech(start:start+(l-overlap))=cleanspeech(i,overlap:end);
    start=start+l-overlap-1;
end
cleanSpeech(start:length(y))=cleanspeech(totseg,1:(length(y)-start)+1);
cleanSpeech=cleanSpeech(1:length(y));
   

end

