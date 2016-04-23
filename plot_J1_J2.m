%plots all J1,J2,nq values for all types of noise and all values of Q
clear all;
close all;
type = ['white ';'babble';'train '];
orderType = 'estimated';
loadFromPath = ['Results\Rnew all noise ',orderType,' order\'];
noise = cellstr(type);
noisename = '';
fig = 1;
toDrawJ1 = zeros(2,10);
toDrawJ2 = zeros(2,10);
toDrawnq = zeros(2,10);
rgb = zeros(1,3);
allPossibleColors = randperm(255);

for i = 1:3
    noisename = char(noise(i));
    load([loadFromPath,noisename,'\','J1_J2_plots_',noisename,'.mat']);
    for j = 1:size(J1v,1)
        figure(fig);
        m = 1;
        l = 1;
        p = 1;
        h = [];
        strLabel = cellstr('');
        for k = 1:size(J1v,2)
            toDrawJ1(1,:) = J1v(j,k,:);
            toDrawJ1(2,:) = J1s(j,k,:);
            toDrawJ2(1,:) = J2v(j,k,:);
            toDrawJ2(2,:) = J2s(j,k,:);
            toDrawnq(1,:) = nqv(j,k,:);
            toDrawnq(2,:) = nqs(j,k,:);
            rgb = allPossibleColors(m:m+2)/255;
            h(l) = plot(toDrawnq(1,:),toDrawJ1(1,:),'-o','Color',rgb);hold on;grid on;
            strLabel{l} = (['J1 voiced Q',num2str(k-1)]);
            set(h(l),'LineWidth',1.25,'MarkerSize',4);
            h(l+1) = plot(toDrawnq(1,:),toDrawJ2(1,:),'-s','Color',rgb);hold on;grid on;
            set(h(l+1),'LineWidth',1.25,'MarkerSize',4);
            strLabel{l+1} = (['J2 voiced Q',num2str(k-1)]);
            
            rgb = allPossibleColors(m+3:m+5)/255;
            h(l+2) = plot(toDrawnq(2,:),toDrawJ1(2,:),'-o','Color',rgb);hold on; grid on;
            set(h(l+2),'LineWidth',1.25,'MarkerSize',4);
            strLabel{l+2} = (['J1 silent Q',num2str(k-1)]);
            h(l+3) = plot(toDrawnq(2,:),toDrawJ2(2,:),'-s','Color',rgb);hold on;grid on;
            set(h(l+3),'LineWidth',1.25,'MarkerSize',4);
            strLabel{l+3} = (['J2 silent Q',num2str(k-1)]);
            m = m+6;
            l = l+4;
        end
        xlabel('n_q');ylabel('J_1, J_2');
        title([num2str((j-1)*5),'dB',' ',noisename]);
        legend(h,strLabel);
        saveas(figure(fig),[loadFromPath,noisename,'\','J1_J2 plot for',num2str((j-1)*5),'dB.fig']);
        saveas(figure(fig),[loadFromPath,noisename,'\','J1_J2 plot for',num2str((j-1)*5),'dB.jpg']);
        fig = fig+1; 
    end
end

    

