function[] = show_fit(filename_body,filename_tail)

load(filename_body,'Xfit','X');
Xfit_all = Xfit; Xall = X;
load(filename_tail,'Xfit','X');
Xfit_all = cat(1,Xfit_all,Xfit); 
Xall = cat(1,Xall,X);
Nframe = size(Xfit_all,3);

figure; 
subplot(1,1,1); 
for n = 1:Nframe
    plot(Xall(:,1,n),Xall(:,2,n),'.k','MarkerSize',14);
    hold on;
    plot(Xfit_all(:,1,n),Xfit_all(:,2,n),'ob','MarkerSize',8,'LineWidth',2);
    xlim([0 1500]); ylim([0 1500]);
    title(num2str(n));
    drawnow; pause(0.025);
    cla;
end

%% run function
%MIKE (change filenames to fitted data for body and tail):
filename_body =  'C:\PhD 2nd Year\DLC Tracking Data\Extinction\mouse10_extinction_p1_2024-10-04-155321-0000DLC_resnet50_Fear Extinction No ImplantOct15shuffle1_500000_body_fit';
filename_tail = 'C:\PhD 2nd Year\DLC Tracking Data\Extinction\mouse10_extinction_p1_2024-10-04-155321-0000DLC_resnet50_Fear Extinction No ImplantOct15shuffle1_500000_tail_fit';
show_fit(filename_body,filename_tail)
