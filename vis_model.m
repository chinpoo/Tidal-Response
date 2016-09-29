% func: draw circles to visualize radial boundaries in the model
% % input args: 
%     r0:         column vector of radius of material layers, from small to large radius
%     r_lower:    column vector of radius of lh lower boundaries
%     r_upper:    column vector of radius of lh upper boundaries

function vis_model(r0,r_lower,r_upper)
    figure;
    for i = 1:length(r0)
        draw_circle(r0(i),1);          % black dashed circle for material layer boundary
        hold on;
    end
    for i = 1:length(r_lower)
        draw_circle(r_lower(i),2);     % red solid circle for lh layer lower boundary
        hold on;
        draw_circle(r_upper(i),3);     % green solid circle for lh layer upper boundary
        hold on;
    end
    axis square;
    hold off;

    % nested func: draw a circle
    function draw_circle(r,ind)
        theta = 0:0.01:2*pi;
        x = r*cos(theta);
        y = r*sin(theta);
        c = {'black','red','green'};
        s = {':','-','-'};
        w = [1,1.5,1.5];
        p = plot(x,y,s{ind});
        set(p,'LineWidth',w(ind),'Color',c{ind});
    end

end
