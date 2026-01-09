function[] = plot3d_video(X, make_lines)
[Np,~,Ntt] = size(X);
figure; 
for t = 1:1:Ntt
    hold on;
    for m = 1:Np
        plot3(squeeze(X(m,1,t)),squeeze(X(m,2,t)),squeeze(X(m,3,t)),'.','MarkerSize',20);
    end
    %
    if make_lines
        line([X(1,1,t) X(3,1,t)],[X(1,2,t) X(3,2,t)],[X(1,3,t) X(3,3,t)],'Color','k','LineWidth',2);
        line([X(2,1,t) X(3,1,t)],[X(2,2,t) X(3,2,t)],[X(2,3,t) X(3,3,t)],'Color','k','LineWidth',2);
        %
        line([X(3,1,t) X(7,1,t)],[X(3,2,t) X(3,2,t)],[X(3,3,t) X(7,3,t)],'Color','k','LineWidth',2);
        line([X(3,1,t) X(8,1,t)],[X(3,2,t) X(8,2,t)],[X(3,3,t) X(8,3,t)],'Color','k','LineWidth',2);
        line([X(7,1,t) X(9,1,t)],[X(7,2,t) X(9,2,t)],[X(7,3,t) X(9,3,t)],'Color','k','LineWidth',2);
        line([X(8,1,t) X(10,1,t)],[X(8,2,t) X(10,2,t)],[X(8,3,t) X(10,3,t)],'Color','k','LineWidth',2);
        %
        line([X(7,1,t) X(8,1,t)],[X(7,2,t) X(8,2,t)],[X(7,3,t) X(8,3,t)],'Color','k','LineWidth',2);
        line([X(9,1,t) X(10,1,t)],[X(9,2,t) X(10,2,t)],[X(9,3,t) X(10,3,t)],'Color','k','LineWidth',2);
        line([X(9,1,t) X(11,1,t)],[X(9,2,t) X(11,2,t)],[X(9,3,t) X(11,3,t)],'Color','k','LineWidth',2);
        line([X(10,1,t) X(12,1,t)],[X(10,2,t) X(12,2,t)],[X(10,3,t) X(12,3,t)],'Color','k','LineWidth',2);
        line([X(11,1,t) X(12,1,t)],[X(11,2,t) X(12,2,t)],[X(11,3,t) X(12,3,t)],'Color','k','LineWidth',2);
        %
        line([X(3,1,t) X(4,1,t)],[X(3,2,t) X(4,2,t)],[X(3,3,t) X(4,3,t)],'Color','k','LineWidth',2);
        line([X(4,1,t) X(5,1,t)],[X(4,2,t) X(5,2,t)],[X(4,3,t) X(5,3,t)],'Color','k','LineWidth',2);
        line([X(5,1,t) X(6,1,t)],[X(5,2,t) X(6,2,t)],[X(5,3,t) X(6,3,t)],'Color','k','LineWidth',2);
    end
    %
    if make_lines
        xlim(0.6*[-15 15]); ylim(0.6*[-15 15]); zlim(0.6*[-15 15]);
    else
        xlim([-25 25]); ylim([-25 25]); zlim([-5 20]);
    end
    view([-70 20]); title(num2str(t));
    % view([45 0]); 
    drawnow; 
    pause(0.04); 
    if t ~= Ntt
        cla;
    end
end