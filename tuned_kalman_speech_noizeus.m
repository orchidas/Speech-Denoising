function [cleanSpeech] = tuned_kalman_speech_noizeus(filename, noiseType, orderType)
%Applies tuned Kalman filter with order estimation on noisy speech.
%filename - name of .wav speech file from NOIZEUS
%noiseType - white, train or babble
%orderType = estimated or fixed

parentpath = fileparts(pwd);
SNR = [0,5,10];

%this folder contains MATLAB files needed to calculate PESQ
%download it from http://ecs.utdallas.edu/loizou/speech/software.htm
%and extract it in the parent directory
addpath(strcat(parentpath,'\composite\'));

%this folder contains clean and corrupted .wav speech files downloaded from
%NOIZEUS database - http://ecs.utdallas.edu/loizou/speech/noizeus/
soundpath = strcat(parentpath,'\Noisy speech samples\');

%folder where results are saved - create if does not exist
saveToPath = ['Results\Rnew all noise ',orderType, ' order\',noiseType,'\',filename,'\'];
if exist(saveToPath, 'dir') == 0
    mkdir(saveToPath);
end
   
%writing results to txt file
[fileID, message] = fopen([saveToPath,filename,'_',noiseType,'_results.txt'],'w+');
fprintf(fileID,'%s %s %s %s %s %s %s %s\r\n','R(new)','Q_chosen','log(Q)',...
    'SNR','SegSNR_before','SegSNR_after','PESQ','Average_order');

for snri = 1:length(SNR)
    
    %read a clean audio signal
    [z,fs] = audioread(strcat(soundpath,'clean\',filename,'.wav'));
    %read a corrupted audio signal
    if(strcmp(noiseType, 'white') == 1)
        wn = audioread([soundpath,'white_noise.wav']);
        [noise,snr] = makeSNR(z, wn, SNR(snri));
        y = noise + z;
    else
        [y,fs] = audioread(strcat(soundpath, noiseType,'\', num2str(SNR(snri)),...
            'dB\',filename,'_',type,'_sn',num2str(SNR(snri)),'.wav'));
    end
    y=y';
    z=z';
    
    %dividing into 80ms frames with 10ms overlap
    start=1;
    l=0.08*fs;
    overlap=0.01*fs;
    totseg=ceil(length(y)/(l-overlap));
    segment=zeros(totseg,l);
    zseg=zeros(totseg,l);
    
    for i=1:totseg-1
        segment(i,1:l)=y(1,start:start+l-1);
        zseg(i,1:l)=z(1,start:start+l-1);
        start=(l-overlap)*i+1;
    end
    
    segment(totseg,1:length(y)-start+1)=y(start:length(y));
    zseg(totseg,1:length(z)-start+1)=z(start:length(z));
    cleanspeech=zeros(totseg,l);
    cleanSpeech=zeros(1,length(y));
    
    %determine order
    if strcmp(orderType,'fixed') == 1
        order = ones(totseg,1).*15;
    else
        order = findOrder(segment,SNR(snri),type,saveToPath);
    end
    
    %calculate measurement noise variance R
    R = measurementNoiseNew(segment, fs);
    
    J1=zeros(1,10);
    J2=zeros(1,10);
    nq=zeros(1,10);
    u=1;
    
    for m = 0:6
        
        segsnr_before=0;
        segsnr_after=0;
        Q_arr=zeros(1,totseg);
        
        for i=1:totseg
            %initializing
            X=y(1:order(i))';
            P=zeros(l,order(i),order(i));
            t1=zeros(order(i),order(i));
            H=[zeros(1,order(i)-1),1];
            G=H';
            
            %first iteration of Kalman filter
            [A,Q1]=lpc(segment(i,:),order(i));
            temp=eye(order(i));
            PHI=[temp(2:order(i),:);-fliplr(A(2:end))];
            
            if(i == 1)
                P(1,:,:)=R*eye(order(i));
            else
                P(1,:,:) = Y(:,:);
            end
            
            %calculating optimum value of process noise variance, Q
            q=1;
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
            
            %interpolate nq, J1 and J2 to increase resolution, and to get more
            %accurate approximation of Q
            nqi = -5:0.25:4;
            J2i = interp1(nq,J2,nqi);
            J1i = interp1(nq,J1,nqi);
            [nq_nom,Jc]=intersections(nqi,J1i,nqi,J2i);
            
            if m < 3
                J2_desired = (0.25*(m+1)*(Jc-min(J2i))) + min(J2i);
            else
                J2_desired = (0.25*(m-3)*(max(J2i)-Jc))+ Jc;
            end
            
            [difference, index] = min(abs(J2i - J2_desired));
            Q = 10^(nqi(index));
            
            %plot J1,J2 for voiced frame
            if(i == 4)
                figure(1);
                plot(nqi,J1i,'-b+');hold on;grid on;
                plot(nqi,J2i,'-b*');hold on;grid on;
                scatter(nqi(index),J2i(index),'k');
            end
            
            %plot J1,J2 for silent frame
            if(i == totseg - 1)
                figure(1);
                plot(nqi,J1i,'-r+');hold on;grid on;
                plot(nqi,J2i,'-r*');hold on;grid on;
                scatter(nqi(index),J2i(index),'k');
                xlabel('nq','FontSize',18);
                ylabel('JI,J2','FontSize',16);
                axis([min(nqi)-2, max(nqi), 0, 1]);
                hold off;
                legend('J_1','J_2');
            end
            
            Q_arr(u)=Q;
            u=u+1;
            
            for j=1:length(segment(i,:))
                X_=PHI*X;
                t1(:,:)=P(j,:,:);
                P_=(PHI*t1*PHI')+(G*Q*G');
                K=(P_*H')*(inv(H*P_*H'+R));
                t1=(eye(order(i))-K*H)*P_;
                P(j+1,:,:)=t1(:,:);
                e=segment(i,j)-(H*X_);
                X=X_+K*e;
                cleanspeech(i,j)=X(end);
                
            end
            
            %adjust a posteriori error covariance matrix dimensions
            if(i< totseg)
                t2 = zeros(order(i),order(i));
                t2(:,:) = P(j-1,:,:);
                Y = adjustDimensions(t2, order(i+1));
            end
            
            %second iteration of Kalman filter with lpc calculated from
            %cleaned speech
            [A,Q]=lpc(cleanspeech(i,:),order(i));
            PHI=[temp(2:order(i),:);-fliplr(A(2:end))];
            X=cleanspeech(i,1:order(i))';
            
            if i==1
                P0=R*eye(order(i));
            else
                P0 = Z(:,:);
            end
            
            for j=1:length(segment(i,:))
                X_=PHI*X;
                P_=(PHI*P0*PHI')+(G*Q*G');
                K=(P_*H')*(inv(H*P_*H'+R));
                P0=(eye(order(i))-K*H)*P_;
                e=segment(i,j)-(H*X_);
                X=X_+K*e;
                cleanspeech(i,j)=X(end);
            end
            
            if(i< totseg)
                Z = adjustDimensions(P0, order(i+1));
            end
            
            %calculate Segmental SNR
            segsnr_before=segsnr_before+log10(rms(zseg(i,:))/rms(zseg(i,:)-segment(i,:)));
            segsnr_after=segsnr_after+log10(rms(zseg(i,:))/rms(zseg(i,:)-cleanspeech(i,:)));    
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
        
        %normalizing
        z=z./abs(1.2*max(z));
        y=y./abs(1.2*max(y));
        cleanSpeech=cleanSpeech./abs(1.2*max(cleanSpeech));
        
        %qualitative measure of noise removed
        figure(2);ylabel('Normalised amplitude');xlabel('Time in seconds');
        subplot(3,1,1);plot((1:length(z))/fs,z);title('original speech');axis([0, length(z)/fs, -1, 1]);
        subplot(3,1,2);plot((1:length(y))/fs,y,'k');title('corrupted speech');axis([0, length(y)/fs, -1, 1]);
        subplot(3,1,3);plot((1:length(cleanSpeech))/fs,cleanSpeech,'r');title('cleaned speech');axis([0, length(cleanSpeech)/fs, -1, 1]);
        soundsc(y,fs);
        soundsc(cleanSpeech,fs);
        audiowrite(cleanSpeech,fs,strcat(saveToPath,filename,'_Q',num2str(m),'_',noiseType,'_sn',num2str(SNR(snri)),'_enhanced.wav'));
        saveas(figure(2),[saveToPath,'Waveform_',filename,'_Q',num2str(m),'_',noiseType,'_sn',num2str(SNR(snri))]);
        
        figure(3);
        subplot(3,1,1);spectrogram(z,64,16,1024,fs,'yaxis');
        title('Spectrum of clean speech');
        subplot(3,1,2);spectrogram(y,64,16,1024,fs,'yaxis');
        title('Spectrum of noisy speech');
        subplot(3,1,3);spectrogram(cleanSpeech,64,16,1024,fs,'yaxis');
        title('Spectrum of cleaned speech after processing');
        saveas(figure(3), [saveToPath,'Spectrogram_',filename,'_Q',num2str(m),'_',noiseType,'_sn',num2str(SNR(snri))]);
        
        %quantitative measure of noise removed - Seg SNR
        disp('The segmental snr before processing is :');
        segsnr_before = 20*segsnr_before/totseg
        disp('The segmental snr after processing is :');
        segsnr_after = 20*segsnr_after/totseg
        
        %PESQ
        psq = pesq(fs,strcat(soundpath,'clean\',filename,'.wav'),...
            strcat(saveToPath, filename,'_Q',num2str(m),'_',noiseType,'_sn',num2str(SNR(snri)),'_enhanced.wav'));
        fprintf(fileID,'%f %s %s %d %f %f %f %d\r\n', R, ['Q',num2str(m),'=',num2str(mean(Q_arr))],...
            ['n=',num2str(log10(mean(Q_arr)))], SNR(snri), segsnr_before, segsnr_after, psq, round(mean(order)));
        close all;  
    end
end
fclose('all');
end


