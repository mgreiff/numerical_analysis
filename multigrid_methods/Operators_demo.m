npts = 7;
xx = 1:npts;
rMat = reshape(wgn(1,npts.^2,1),npts,npts)/10;
pMat = FMGprolong(rMat);

if 1
    subplot(1,2,1)
    hold on;
    
    [X,Y]=meshgrid(linspace(0,1,length(pMat)));
    s1 = surf(X,Y,pMat,'EdgeColor','r')
    alpha(s1,.5)
    shading interp;
    
    surf(X,Y,pMat,'FaceColor','None','EdgeColor','r')
    plot3(X,Y,pMat,'ro')
    
    [X,Y]=meshgrid(linspace(0,1,length(rMat)));
    surf(X,Y,rMat,'FaceColor','None','EdgeColor','k','LineWidth',2)
    
    title('The prolong operator.')
    axis([0,1,0,1,-1,1])
    xlabel('x')
    ylabel('y')
    zlabel('u(x,y)')
    hold off;
end

sMat = FMGrestrict(rMat);

if 1
    subplot(1,2,2)
    hold on;
    [X,Y]=meshgrid(linspace(0,1,length(sMat)));
    s1 = surf(X,Y,sMat,'EdgeColor','b', 'LineWidth',2)
    alpha(s1,.5)
    shading interp;

    surf(X,Y,sMat,'FaceColor','None','EdgeColor','b')
    plot3(X,Y,sMat,'bo')
    
    [X,Y]=meshgrid(linspace(0,1,length(rMat)));
    surf(X,Y,rMat,'FaceColor','None','EdgeColor','k','LineWidth',2)
    axis([0,1,0,1,-1,1])
    title('The restric operator.')
    xlabel('x')
    ylabel('y')
    zlabel('u(x,y)')
    hold off;
end
