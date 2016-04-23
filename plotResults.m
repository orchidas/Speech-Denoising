%this script plots segSNR and PESQ values for each type of noise, of fixed
%or variable order.
clf;
clear all;
otype = ['fixed    ';'estimated'];
orderType = cellstr(otype);
ntype = ['babble';'train ';'white '];
noiseType = cellstr(ntype);
%fileparts(pwd) gets parent directory of current directory
saveToPath = [fileparts(pwd),'\Report\figures\'];
fig  = 1;

for i = 1:2
    order = char(orderType(i));
    for j = 1:3
        noise = char(noiseType(j));
        readFromPath = ['Results\Rnew all noise ',order,' order\',noise,'\'];
        filename = [readFromPath,'sp12_',noise,'_results.txt'];
        fileID = fopen(filename,'r');
        formatSpec = '%s';
        N = 8;
        C_text = textscan(fileID,formatSpec,N,'Delimiter',' ');
        C_data = textscan(fileID,'%f %s %s %d %f %f %f %d');
        Qstr = char(C_data{2}{:});
        Qstr = strtrim(Qstr(:,4:end));
        nstr = char(C_data{3}{:});
        nstr = strtrim(nstr(:,3:end));
        Q = zeros(1,size(Qstr,1));
        n = zeros(1,size(nstr,1));
        for k = 1:length(Q)
            Q(k) = str2double(Qstr(k,:));
            n(k) = str2double(nstr(k,:));
        end
        segSNR = C_data{6};
        PESQ = C_data{7};
        
        %plot segSNR
        figure(fig);
        h = plot(n(1:7),segSNR(1:7),'b--o',n(8:14), segSNR(8:14),'g--o', n(15:21), segSNR(15:21),'r--o');
        xlabel('log_{10}(Q)','FontSize',14);
        ylabel('Seg SNR (dB)','FontSize',14);
        set(h,'LineWidth',2,'MarkerSize',6);
        grid on;
        h_leg = legend('SNR 0dB','SNR 5dB','SNR 10dB');
        set(h_leg,'FontSize',12);
        set(gca,'FontSize',12);
        print([saveToPath,'segSNR_',order,'_',noise],'-depsc');
        
        %plot PESQ
        figure(fig+1);
        h = plot(n(1:7),PESQ(1:7),'b--o',n(8:14), PESQ(8:14),'g--o', n(15:21), PESQ(15:21),'r--o');
        xlabel('log_{10}(Q)','FontSize',14);
        ylabel('PESQ','FontSize',14);
        set(h,'LineWidth',2,'MarkerSize',6);
        grid on;
        h_leg = legend('SNR 0dB','SNR 5dB','SNR 10dB');
        set(h_leg,'FontSize',12);
        set(gca, 'FontSize',12);  
        print([saveToPath,'pesq_',order,'_',noise],'-depsc');
        
        fig = fig+2;
    end
end



