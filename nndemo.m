function nndemo

x=rand(1,1000);
y=rand(1,1000);
nntest(x,y,1,1)

x=rand(1,1000);
y=x; y(2:end)=x(1:end-1)*0.5;
nntest(x,y,1,1)