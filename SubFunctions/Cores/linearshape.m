function [N,dN]=linearshape(dx,Lx)
    c1  = (abs(dx)<=Lx);
    N1 = 1-abs(dx)/Lx; 
    dN1 = -sign(dx)/Lx;
    N = N1.*c1; 
    dN = c1.*dN1;
