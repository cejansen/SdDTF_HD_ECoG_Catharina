function dotplot(IN,OUT,TOT,res,scale)
%IN=INgs;OUT=OUTgs;TOT=TOTgs;
[X,Y] = meshgrid(1:9);
figure; hold on
plot(X,Y,'k');
plot(Y,X,'k'); axis off
n = 100;N=1:63;
cmap = parula(n);colorbar('Ticks',[0 1],'TickLabels',{'More outgoing','More Ingoing'})
a=[1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8];
b=[1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7];
mo=max(OUT);mi=max(IN);
if max(mo>mi)
    nOUT=OUT./mo;
    nIN=IN./mo;
else
    nOUT=OUT./mi;
    nIN=IN./mi;
end
ninout=((nIN-abs(nOUT))+1)/2;
for i = 1:63
    color=cmap(round(100.*ninout(i)),:);
    circle2(a(i)+0.5,b(i)+0.5,TOT(i)*scale,[color 0.7]);
    rescheck=any(res(:) == i);
    if rescheck==1
    text(a(i)+0.5,b(i)+0.5,num2str(N(i)),"FontWeight","bold");
    else
    text(a(i)+0.5,b(i)+0.5,num2str(N(i)));
    end
end
end

function h = circle2(x,y,r,color)
d = r*2;
px = x-r;
py = y-r;
h = rectangle('Position',[px py d d],'Curvature',[1,1],'FaceColor',color,'LineStyle','none');
daspect([1,1,1])
end
