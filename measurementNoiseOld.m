function [R] = measurementNoiseOld(x,fs)
%old method of determining R according to spectral energy as 
%given in CMI paper

totseg = size(x,1);
or_num = zeros(1,totseg);
or_den = zeros(1,totseg);
ratio_or = zeros(1,totseg);
count1 = 0;
count2 = 0;

for i=1:totseg
     
   [S0,f0,t0,p0]=spectrogram(x(i,:),64,16,1024,fs);
   s=size(S0);
 
       for j=1:s(1)
       %for components below 2kHz(speech)
        if(f0(j)<=2000)
            for k=1:s(2)
            or_num(i)=or_num(i)+((abs(S0(j,k))).^2);
            end
            count1=count1+1;
        %for components above 2kHz(noise)
        else
            for k=1:s(2)    
            or_den(i)=or_den(i)+((abs(S0(j,k))).^2);
            end
            count2=count2+1;
        end
       end
       or_num(i)=or_num(i)/(count1*s(2));
       or_den(i)=or_den(i)/(count2*s(2));
       ratio_or(i)=or_num(i)/or_den(i);
       count1=0;
       count2=0;
       if(ratio_or(i)<=1||i==1||i==2||i==totseg||i==totseg-1)
            u(k)=var(x(i,:));
        k=k+1;
       else
            continue;
       end
end
R=sum(u)/k


end

