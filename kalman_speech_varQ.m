function [cleanSpeech] = kalman_speech_varQ(x,fs)

%Tuned Kalman filter for speech enhancement with different values of
%process noise covariance

%Inputs:
%x - noisy speech sample
%fs - sampling rate
%Output:
%cleanSpeech - enhanced speech

prompt='What is the value of Q? 0 for Qnom(from lpc),1 for Q1,2 for Q2,3 for Qc,4 for Q3,5 for Q4      ';
in=input(prompt);

while in<0 || in>5
    disp('Wrong input,input must be between 0 and 5');
    in=input(prompt);
    if in>=0 && in<=5
        break;
    end
end
    
Q_pos=[-3,-0.7,0,1,3];

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


R = measurementNoiseNew(segment,fs);
J1=zeros(1,10);
J2=zeros(1,10);
nq=zeros(1,10);
u=1;

X=y(1:p)';
P=zeros(l,p,p);
P(1,:,:)=R*eye(p);
t1=zeros(p,p);
Q_arr=zeros(1,totseg);

for i=1:totseg
    
    %first iteration of Kalman filter
    [A,Q1]=lpc(segment(i,:),p);
    temp=eye(p);
    PHI=[temp(2:p,:);-fliplr(A(2:end))];
    
    %tuning the filter by calculating optimum value of process noise
    %variance
   q=1;
   if i~=1 && in~=0

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
     
     
     [nq_nom,J]=intersections(nq,J1,nq,J2);
     
     %sensitivity metrics gets higher preference and we select Q<Qcomp
     if numel(nq_nom)~=0
            Q=10^(nq_nom+Q_pos(in));
     else
       Q=Q1;
     end
     
   else
       Q=Q1;
   end
   Q_arr(u)=Q;
   u=u+1;
   
   for j=1:length(segment(i,:))
        X_=PHI*X;
        t1(:,:)=P(j,:,:);
        P_=(PHI*t1*PHI')+(G*Q*G');
        K=(P_*H')*(inv(H*P_*H'+R));
        t1=(eye(p)-K*H)*P_;
        P(j+1,:,:)=t1(:,:);
        e=segment(i,j)-(H*X_);
        X=X_+K*e;
        cleanspeech(i,j)=X(end);
        
   end
    P(1,:,:)=P(j-1,:,:);
       
    %second iteration of Kalman filter with lpc calculated from
    %cleaned speech
    
    [A,Q]=lpc(cleanspeech(i,:),p);
    PHI=[temp(2:p,:);-fliplr(A(2:end))];
     if i==1
        X=cleanspeech(i,1:p)';
        P0=temp*R;
     end
    
    for j=1:length(segment(i,:))
        
        X_=PHI*X;
        P_=(PHI*P0*PHI')+(G*Q*G');
        K=(P_*H')*(inv(H*P_*H'+R));
        P0=(eye(p)-K*H)*P_;
        e=segment(i,j)-(H*X_);
        X=X_+K*e;
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

