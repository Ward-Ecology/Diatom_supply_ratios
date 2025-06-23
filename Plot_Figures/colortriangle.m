function colortriangle(n)
%%

c=0;
for i=0:n
    for j=0:(n-i)
        c=c+1;
        
        x(c) = i/n + j/(2*n);
        y(c)  = j/n * sqrt(3/4);

        r(c) = max(min(i./n,1),0);
        g(c) = max(min(j./n,1),0);
        b(c) = max(min(1 - r(c) - g(c),1),0);
    end
end
%

[vMGFPoints, vMGFcells] = voronoin([x' y']);
[vMGFx, vMGFy] = voronoi(x,y);


hold on
rgb = [r' g' b'].^0.5;

rgb=rgb./max(rgb,[],2).*0.8;
for i = 1:c
    xx = vMGFPoints(vMGFcells{i},1);
    yy = vMGFPoints(vMGFcells{i},2);

    patch(xx,yy, rgb(i,:),'EdgeColor','k','LineW',1)
end
for i = 1:c
    patch(vMGFPoints(vMGFcells{i},1),vMGFPoints(vMGFcells{i},2), rgb(i,:),'EdgeColor','none')
end
axis equal


axis off

text(0,0,'N','HorizontalAlignment','right')
text(1,0,'Fe','HorizontalAlignment','left')
text(0.5,sqrt(3/4),'Si','HorizontalAlignment','center','VerticalAlignment','bottom')
text(0.25,0.5.*sqrt(3/4),'N,Si ','HorizontalAlignment','right','VerticalAlignment','middle')
text(0.75,0.5.*sqrt(3/4),' Si,Fe','HorizontalAlignment','left','VerticalAlignment','middle')
text(0.5,00,'N,Fe','HorizontalAlignment','center','VerticalAlignment','top')


plot([0 1 0.5 0],[0 0 sqrt(3/4) 0],'k-','LineW',1)

%%
end