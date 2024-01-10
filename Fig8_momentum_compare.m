% comparing between
% 0% tau 100% heat, 100% tau x 100% heat, 100% tau x 0% heat, 0% tau x 0%
% heat

bye; clf
load Yixi_c_PrettyFive.mat
Yixi_c_PrettyFive(2,:)=[];
day_start=[51.8297323193263];

%% loading 
n=1;
% OutputFile='PWP_output_1_1_0_0tau';
OutputFile='PWP_output_ControlFluxes_NoMomentum';
load(OutputFile)
ml(n,:)=pwp_output.ml;

n=n+1;
% OutputFile='PWP_output_1_1_1_0tau';
OutputFile='PWP_output_NoFlux_NoMomentum';
load(OutputFile)
ml(n,:)=pwp_output.ml;

n=n+1;
OutputFile='PWP_output_ControlFluxes';
load(OutputFile)
ml(n,:)=pwp_output.ml;

n=n+1;
OutputFile='PWP_output_NoFlux';
load(OutputFile)
ml(n,:)=pwp_output.ml;

%% plotting
pwp_output.day=pwp_output.time+day_start;

set(gca,'FontSize',15);
set(gcf,'PaperPositionMode','auto')

ylabel('Mixed-layer depth (m)','fontsize',17);

xlabel('Year Day','fontsize',23)
hold on

for i=1:size(ml,1)
    plot(pwp_output.day,ml(i,:),'color',Yixi_c_PrettyFive(i,:),'LineWidth',2)
end


xticks([50:5:110]);

% grid minor; box on;
xlim([min(pwp_output.day) max([pwp_output.day])])
ylim([0 247])
grid on; set(gca,'layer','top');
grid minor; set(gca,'layer','top')
box on; set(gca,'layer','top')
ylim([0 225])