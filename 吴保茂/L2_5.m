x=unifrnd(0,12,[1,10000000]); %以0和12为上下限，产生尺寸为1，1000000的矩阵
y=unifrnd(0,9,[1,10000000]);
pinshu=sum(y<x.^2 & x <=3)+sum(y<12-x & x>=3);
area_appr= 12*9*pinshu/10^7