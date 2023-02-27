function h=plotDoubletResponse(SYSOL,SYSCL,TIMEEND)

no_output=size(SYSOL,1);
no_input=size(SYSOL,2);

% Setup simulation
t=0:.01:TIMEEND;
doublet=(t>1)-2*(t>2)+(t>3);

yol=zeros(length(t),no_output,no_input);
ycl=zeros(length(t),no_output,no_input);
for i=1:no_input
    [yol(:,:,i),~]=lsim(SYSOL(:,i),doublet,t);
    [ycl(:,:,i),~]=lsim(SYSCL(:,i),doublet,t);
end


% Plot
plot_colors = [55, 126, 184; ...
              228,  26,  28; ...
               77, 175,  74; ...
              152,  78, 163; ...
              255, 127,   0; ...
              255, 255,  51]/255;
h=figure;
no_fig=1;
for i=1:no_output
    for j=1:no_input
        subplot(no_output,no_input,no_fig); hold on
        plot(t,yol(:,i,j),'LineWidth',1.5,'Color',plot_colors(1,:));
        plot(t,ycl(:,i,j),'LineWidth',1.5,'Color',plot_colors(2,:));
        if i==1
            title(SYSOL.InputName{j});
            legend('OL','CL', ...
                   'Location','southeast', ...
                   'NumColumns',2);
        end
        ylabel(append('$',SYSOL.OutputName{i},'$'),'Interpreter','latex','FontSize',14);
        if i==no_output
            xlabel('t','Interpreter','latex','FontSize',14);
        end
        hold off
        grid on
        no_fig=no_fig+1;
    end
end
end