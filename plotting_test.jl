using PlotlyJS

x=1:3
y=1:3

p=plot(x,y,Layout(xaxis=attr(range=[-2,2],maxallowed=10)))
savefig(p,"test.html")
