function plot_lapMask(lapmask, chanlocs)

    X = [chanlocs.X];
    Y = [chanlocs.Y];
    Z = [chanlocs.Z];
    
    figure()
    scatter3(Y,X,Z);
    hold on;
    text(Y,X,Z,{chanlocs.labels})
    
    pbaspect([1 1 1]);
    [cs,cd] = find(lapmask-eye(length(lapmask))~=0);
    for i = 1:length(cs)
        plot3([Y(cs(i)),Y(cd(i))],[X(cs(i)),X(cd(i))],[Z(cs(i)),Z(cd(i))]);
    end
