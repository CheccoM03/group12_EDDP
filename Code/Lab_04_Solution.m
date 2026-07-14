clear all
close all
clc

set(0,'defaultaxesfontsize',16);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultStemLineWidth', 1.5);
%%

[n,p]=uigetfile('.mat','choose the data file');
cd(p);
load(n);

% Plot the data
dt=1/fsamp;
N=length(Data);
t=0:dt:(N-1)*dt;

figure (1)
ax1=subplot(2,1,1);
plot(t,Data(:,1));
set(gca,'fontsize',14)
xlabel('t [s]');
ylabel('Force [N]')
grid on
axis tight

ax2=subplot(2,1,2);
plot(t,Data(:,2),t,Data(:,3));
set(gca,'fontsize',14)
xlabel('t [s]');
ylabel('Acc [m/s^2]')
legend('Accelerometer 17','Accelerometer 18')
axis tight
grid on

linkaxes([ax1 ax2],'x');

%% Calculating of H=Y/X for whole time history
[Xt, freq]=fft_n(hanning(N).*Data(:,1),fsamp); % Since I want to plot the spectrum only for the positive frequencies, I define a function ("fft_n", see the end of the script) to avoid repeating the same operations many times
[Y1t, freq]=fft_n(hanning(N).*Data(:,2),fsamp);
[Y2t, freq]=fft_n(hanning(N).*Data(:,3),fsamp);

Ht1=Y1t./Xt;
Ht2=Y2t./Xt;

figure
ax1=subplot(2,1,1);
semilogy(freq,abs(Ht1))
set(gca,'fontsize',14)
title('Accelerometer 17: H=Y/X')
xlabel('Frequenct [Hz]')
ylabel('|H| [m/(s^2*N)]')
axis tight
grid on
ax2=subplot(2,1,2);
plot(freq,angle(Ht1))
set(gca,'fontsize',14)
xlabel('Frequency [Hz]')
ylabel('\phi [rad]')
axis tight
grid on
linkaxes([ax1 ax2],'x');

figure(9)
ax1=subplot(2,1,1);
semilogy(freq,abs(Ht2))
set(gca,'fontsize',14)
title('Accelerometer 18: H=Y/X')
xlabel('Frequency [Hz]')
ylabel('|H| [m/(s^2*N)]')
axis tight
grid on
ax2=subplot(2,1,2);
plot(freq,angle(Ht2))
set(gca,'fontsize',14)
xlabel('Frequency [Hz]')
ylabel('\phi [rad]')
axis tight
grid on
linkaxes([ax1 ax2],'x');

%% Power and cross-spectra

for Time=[30 60]
    SecPoints=round(Time*fsamp);
    N_OL=floor(0.66*(SecPoints));

    Win=hanning(SecPoints);
    [Gxx,frequency]=autocross(Data(:,1),Data(:,1),fsamp,SecPoints,N_OL,Win); % I use a function ("autocross", see the end of the code)  that I defined once for all to avoid repeating the same operations many times



    [Gyy1]=autocross(Data(:,2),Data(:,2),fsamp,SecPoints,N_OL,Win);
    [Gyy2]=autocross(Data(:,3),Data(:,3),fsamp,SecPoints,N_OL,Win);
    [Gxy1]=autocross(Data(:,1),Data(:,2),fsamp,SecPoints,N_OL,Win);
    [Gxy2]=autocross(Data(:,1),Data(:,3),fsamp,SecPoints,N_OL,Win);
    [Gy1x]=autocross(Data(:,2),Data(:,1),fsamp,SecPoints,N_OL,Win);
    [Gy2x]=autocross(Data(:,3),Data(:,1),fsamp,SecPoints,N_OL,Win);

    coherence1=abs(Gxy1).^2./(Gxx.*Gyy1);
    coherence2=abs(Gxy2).^2./(Gxx.*Gyy2);

    %% Power-spectra of the input and the output

    figure (2)
    semilogy(frequency,Gxx)
    hold on

    if Time==60
        legend ('30 s','60 s')
        xlabel('Frequency [Hz]')
        ylabel('Gxx [N^2]')
        title(['Power spectrum of the input'])
        grid on
        axis tight
    end


    figure (3)
    ax1=subplot(2,1,1);
    semilogy(frequency,Gyy1)
    hold on
    if Time==60
        legend ('30 s','60 s')
        xlabel('Frequency [Hz]')
        ylabel('Gyy [(m/s^2)^2')
        title(['Output Power spectrum: point 17'])
        grid on
        axis tight
    end


    ax2=subplot(2,1,2);
    semilogy(frequency,Gyy2)
    hold on
    set(gca,'fontsize',14)
    xlabel('Frequency [Hz]')
    ylabel('Gyy [(m/s^2)^2')
    title(['Output Power spectrum: point 18'])
    if Time==60
        legend ('30 s','60 s')
    end
    grid on
    axis tight
    linkaxes([ax1 ax2],'x');

    %% Cross-spectra
    figure (4)
    ax1=subplot(2,1,1);
    semilogy(frequency,abs(Gxy1))
    hold on
    if Time==60
        legend ('30 s','60 s')
        xlabel('Frequency [Hz]')
        ylabel('|Gxy| [(m/s^2)N]')
        title(['Cross-spectrum: point 17 - Magnitude'])
        grid on
        axis tight
    end

    ax2=subplot(2,1,2);
    plot(frequency,angle(Gxy1))
    hold on
    if Time==60
        legend ('30 s','60 s')
        xlabel('Frequency [Hz]')
        ylabel('\angle(Gxy) [rad]')
        title(['Cross-spectrum: point 17 - Phase'])
        grid on
        axis tight
    end
    linkaxes([ax1 ax2],'x');

    figure (5)
    ax1=subplot(2,1,1);
    semilogy(frequency,abs(Gxy2))
    hold on
    if Time==60
        xlabel('Frequency [Hz]')
        ylabel('|Gxy| [(m/s^2)N]')
        title(['Cross-spectrum: point 18 - Magnitude'])
        legend ('30 s','60 s')
        grid on
        axis tight
    end

    ax2=subplot(2,1,2);
    plot(frequency,angle(Gxy2))
    hold on
    if Time==60
        legend ('30 s','60 s')
        xlabel('Frequency [Hz]')
        ylabel('\angle(Gxy) [rad]')
        title(['Cross-spectrum: point 18 - Phase'])
        grid on
        axis tight
    end
    linkaxes([ax1 ax2],'x');


    %% Coherence

    figure(6)

    plot(frequency,coherence1);
    hold on

    if Time==60
        legend ('30 s','60 s')
        ylabel('\gamma^2_x_y');
        xlabel('f [Hz]')
        title('Coherence - Point 17')
        grid on
        axis tight
    end




    figure (7)
    plot(frequency,coherence2);
    hold on

    if Time==60
        legend ('30 s','60 s')
        ylabel('\gamma^2_x_y');
        xlabel('f [Hz]')
        title('Coherence - Point 18')
        grid on
        axis tight
    end



    %% H1 ed H2
    H11=Gxy1./Gxx;
    H12=Gxy2./Gxx;

    H21=Gyy1./Gy1x;
    H22=Gyy2./Gy2x;


    if Time==30
        figure (10)
        ax1=subplot(2,1,1);
        semilogy(frequency,abs(H11))
        hold on
        semilogy(frequency,abs(H21))
        ylabel('|H| [m/(s^2*N)]')
        yyaxis right
        plot(frequency,coherence1,'k');
        ylabel('\gamma^2_x_y')
        set(gca,'fontsize',14)
        set(gca,'YColor','Black')
        title(['Accelerometer 17: H_1 and H_2 - ' num2str(Time) 's'])
        xlabel('Frequency [Hz]')
        axis tight
        grid on
        legend('H1','H2','\gamma^2_x_y')
        ax2=subplot(2,1,2);
        plot(frequency,angle(H11))
        hold on
        plot(frequency,angle(H21))
        set(gca,'fontsize',14)
        xlabel('Frequency [Hz]')
        ylabel('\phi [rad]')
        axis tight
        grid on
        legend('H1','H2')
        linkaxes([ax1 ax2],'x');

        figure (11)
        ax1=subplot(2,1,1);
        semilogy(frequency,abs(H12))
        hold on
        semilogy(frequency,abs(H22))
        ylabel('|H| [m/(s^2*N)]')
        yyaxis right
        plot(frequency,coherence2,'k');
        ylabel('\gamma^2_x_y')
        set(gca,'fontsize',14)
        set(gca,'YColor','Black')
        title(['Accelerometer 18: H_1 and H_2 - ' num2str(Time) 's'])
        xlabel('Frequency [Hz]')
        axis tight
        grid on
        legend('H1','H2','\gamma^2_x_y')
        ax2=subplot(2,1,2);
        hold on
        plot(frequency,angle(H12))
        hold on
        plot(frequency,angle(H22))
        set(gca,'fontsize',14)
        xlabel('Frequency [Hz]')
        ylabel('\phi [rad]')
        axis tight
        grid on
        legend('H1','H2')
        linkaxes([ax1 ax2],'x');

        figure (14)
        semilogy(frequency,abs(H11),'b')
        hold on
        semilogy(frequency,abs(H21),'r')


        figure (15)
        plot(frequency,coherence1,'k');
        hold on



        figure (16)
        plot(frequency,angle(H11))
        hold on
        plot(frequency,angle(H21))


        figure (17)
        semilogy(frequency,abs(H12),'b')
        hold on
        semilogy(frequency,abs(H22),'r')


        figure (18)
        plot(frequency,coherence2,'k');
        hold on



        figure (19)
        plot(frequency,angle(H12))
        hold on
        plot(frequency,angle(H22))
    else
        figure (12)
        ax1=subplot(2,1,1);
        semilogy(frequency,abs(H11))
        hold on
        semilogy(frequency,abs(H21))
        ylabel('|H| [m/(s^2*N)]')
        yyaxis right
        plot(frequency,coherence1,'k');
        ylabel('\gamma^2_x_y')
        set(gca,'fontsize',14)
        set(gca,'YColor','Black')
        title(['Accelerometer 17: H_1 and H_2 - ' num2str(Time) 's'])
        xlabel('Frequency [Hz]')
        axis tight
        grid on
        legend('H1','H2','\gamma^2_x_y')
        ax2=subplot(2,1,2);
        plot(frequency,angle(H11))
        hold on
        plot(frequency,angle(H21))
        set(gca,'fontsize',14)
        xlabel('Frequency [Hz]')
        ylabel('\phi [rad]')
        axis tight
        grid on
        legend('H1','H2')
        linkaxes([ax1 ax2],'x');

        figure (13)
        ax1=subplot(2,1,1);
        semilogy(frequency,abs(H12))
        hold on
        semilogy(frequency,abs(H22))
        ylabel('|H| [m/(s^2*N)]')
        yyaxis right
        plot(frequency,coherence2,'k');
        ylabel('\gamma^2_x_y')
        set(gca,'fontsize',14)
        set(gca,'YColor','Black')
        title(['Accelerometer 18: H_1 and H_2 - ' num2str(Time) 's'])
        xlabel('Frequency [Hz]')
        axis tight
        grid on
        legend('H1','H2','\gamma^2_x_y')
        ax2=subplot(2,1,2);
        hold on
        plot(frequency,angle(H12))
        hold on
        plot(frequency,angle(H22))
        set(gca,'fontsize',14)
        xlabel('Frequency [Hz]')
        ylabel('\phi [rad]')
        axis tight
        grid on
        legend('H1','H2')
        linkaxes([ax1 ax2],'x');



        figure (14)
        semilogy(frequency,abs(H11),'b--')
        hold on
        semilogy(frequency,abs(H21),'r--')
        ylabel('|H| [m/(s^2*N)]')


        figure (15)
        plot(frequency,coherence1,'k--')
        hold on
        ylabel('\gamma^2_x_y')
        yyaxis left


        figure (16)
        plot(frequency,angle(H11),'b--')
        hold on
        plot(frequency,angle(H21),'r--')

        figure (17)
        semilogy(frequency,abs(H12),'b--')
        hold on
        semilogy(frequency,abs(H22),'r--')
        ylabel('|H| [m/(s^2*N)]')


        figure (18)
        plot(frequency,coherence2,'k--')
        hold on
        ylabel('\gamma^2_x_y')
        yyaxis left


        figure (19)
        plot(frequency,angle(H12),'b--')
        hold on
        plot(frequency,angle(H22),'r--')
    end

end


figure (14)
xlim([5 15])
set(gca,'fontsize',14)
ylabel('|H| [m/(s^2*N)]')
xlabel('Frequency [Hz]')
grid on
title(['Accelerometer 17: H_1 and H_2 - 30 s and 60 s'])
legend('H1 (30 s)','H2 (30 s)','H1 (60 s)','H2 (60 s)')

figure (15)
xlim([5 15])
set(gca,'fontsize',14)
ylabel('\gamma^2_x_y')
xlabel('Frequency [Hz]')
grid on
legend('\gamma^2_x_y (30 s)','\gamma^2_x_y (60 s)')
title(['Accelerometer 17: \gamma^2_x_y - 30 s and 60 s'])


figure (16)
xlim([5 15])
set(gca,'fontsize',14)
ylabel('\phi [rad]')
xlabel('Frequency [Hz]')
grid on
title(['Accelerometer 17: H_1 and H_2 (\phi) - 30 s and 60 s'])
legend('H1 (30 s)','H2 (30 s)','H1 (60 s)','H2 (60 s)')

figure (17)
xlim([5 20])
set(gca,'fontsize',14)
ylabel('|H| [m/(s^2*N)]')
xlabel('Frequency [Hz]')
grid on
title(['Accelerometer 18: H_1 and H_2 - 30 s and 60 s'])
legend('H1 (30 s)','H2 (30 s)','H1 (60 s)','H2 (60 s)')

figure (18)
xlim([5 20])
set(gca,'fontsize',14)
ylabel('\gamma^2_x_y')
xlabel('Frequency [Hz]')
grid on
legend('\gamma^2_x_y (30 s)','\gamma^2_x_y (60 s)')
title(['Accelerometer 18: \gamma^2_x_y - 30 s and 60 s'])


figure (19)
xlim([5 20])
set(gca,'fontsize',14)
ylabel('\phi [rad]')
xlabel('Frequency [Hz]')
grid on
title(['Accelerometer 18: H_1 and H_2 (\phi) - 30 s and 60 s'])
legend('H1 (30 s)','H2 (30 s)','H1 (60 s)','H2 (60 s)')


%5 Putting the figures in order for the comment section of the laboratory
for ii=19:-1:1
    figure(ii)
end


%% Defined functions
%% Auto / cross average spectrum
%
% [Output, mediacomp, frequencies] = autocross (data1, data2, fsamp,N_win,N_OL,Win)
%
% data1, data2: time histories input
% fsamp: sampling frequency
% N_win: number of points in each subrecords that it has been used to divide the time history
% N_OL: number of points of overlap between two consequent subrecords
% Win: time window used to weight the data


%%
%% This function does the normalisation of the output of the fft function
% Input:
%   - data: input data matrix with size [r,c], with r samples and c signals
%   - fsamp: sampling frequency
% Output:
%   - norm_sp: normalized spectrum (positive frequencies)
%   - freq_vec: frequency vector

function [norm_sp, freq_vec]=fft_n(data,fsamp)

dim=size(data);

if dim(2)>dim(1)
    data=data';
end

N=length(data);
df=fsamp/N;

if (N/2)==(floor(N/2))

    freq_vec=[0:df:(N/2*df)]';
    NF=length(freq_vec);
    sp=fft(data,[],1);
    norm_sp(1,:)=sp(1,:)/N;
    norm_sp(2:N/2,:)=sp(2:N/2,:)/(N/2);
    norm_sp(N/2+1,:)=sp(N/2+1,:)/N;

else

    freq_vec=[0:df:((N-1)/2)*df]';
    NF=length(freq_vec);
    sp=fft(data,[],1);
    norm_sp(1,:)=sp(1,:)/N;
    norm_sp(2:(N+1)/2,:)=sp(2:(N+1)/2,:)/(N/2);

end

end