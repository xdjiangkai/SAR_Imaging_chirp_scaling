clc;clear;
x=1:3;
y=1:5;
temps=[82    79    84
    81    63    84
    80    61    82
    82    65    85
    84    81    86]
[x,y]=meshgrid(x,y);
figure(1);
mesh(x,y,temps);
xlabel('x');
ylabel('y');

figure(2);
xi=1:0.2:3;
yi=1:0.2:5;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,temps,xi,yi,'cubic');
mesh(xi,yi,zi);
xlabel('x');
ylabel('y');
