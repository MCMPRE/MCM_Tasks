c=[1,2,3,4];
c=[c,c]';
a=[1,-1,-1,1;1,-1,1,-1;1,-1,-2,3];
a=[a,-a];
b=[-2,-1,-1/2];
[y,z]=linprog(c,a,b,[],[],zeros(8,1))