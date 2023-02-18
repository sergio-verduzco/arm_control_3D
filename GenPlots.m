% script to plot the results of simulations for the Cb paper.

clear all;
close all;

% FIRST FIGURE
figure;
load('err11d2T8r1.mat');
plot(err11d2T8r1,'r','Linewidth',4);
hold on;
load('err11d2T8r4.mat');
plot(err11d2T8r4,'g','Linewidth',4);
load('err11d2T8r8.mat');
plot(err11d2T8r8,'b','Linewidth',4);
% load('err11d2T8r12.mat');
% plot(err11d2T8r12,'k','Linewidth',4);
% load('err11d2T8r20.mat');
% plot(err11d2T8r20,'c','Linewidth',4);
set(gca,'FontSize',32);
set(gca,'xTick',1:5);
set(gca,'xTickLabel',{'1';'2';'3';'4';'5'});
title('Error (distance to target)','Fontsize',32)
xlabel('time [s]','Fontsize',32,'FontName','Arial')
ylabel('distance [cm]','Fontsize',32,'FontName','Arial')
H = gca;
[LEGH,TEXTH,OUTH,OUTM] = legend(H,'1','4','8');
set(LEGH,'Fontsize',28);

% SECOND FIGURE
clear all;
stray = ['1','2','3','4','5','6','7','8','9','10'];
integrals = zeros(1,8);

for i=1:8
    errs = load(strcat('err11d2T8r',num2str(i),'.mat'));
    name = fieldnames(errs);
    instr = strcat('error = getfield(errs,''',name,''');');
    eval(instr{1});
    intervsD = error.Data(2:end) - error.Data(1:end-1);
    midData = error.Data(1:end-1) + intervsD/2;
    intervsT = error.Time(2:end)-error.Time(1:end-1);
    integrals(i) = sum(midData.*intervsT);
end

figure;
bar(integrals);
set(gca,'FontSize',32);
xlabel('reach attempt','Fontsize',32,'FontName','Arial');
ylabel('error integral [cm s]','Fontsize',32,'FontName','Arial');