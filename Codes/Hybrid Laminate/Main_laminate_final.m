clear
clc
format short g
%MOC_3D_Hybrid
%'7'   - Aligned Nano Hybrid (Aligned Nanofiber + Particulates)
%'9'   - Hybrid Aligned Nanocomposites (Long Fiber + Aligned Nano Fiber)
%'10'  - Hybrid Random Nanocomposites (Long Fiber + Random Nano Fiber)
%'11'  - Hybrid Transversely Aligned Nanocomposites (Long Fiber + Transversely Aligned Nano Fiber)
%'12'  - Hybrid (long fiber + particualates)
a(1)=0;
b(1)=0;

for w=1:19 %increament for angle
    
ort=[ a(w) 0 b(w)] %oreintn of each layer
if (w<3)
a(w+1)=a(w)+5;
b(w+1)=b(w)+5;
else
a=a;
b=b;
end
n_l=3; %number of layers 
order=[10 10 10]; %order of - type of lamina 
layerthick=[.005 .005 .005]; %thickness of each lamina in mm


for(i=1:n_l)
     fibertype=order(i);

     
      [E]=MOC_3D_Hybrid(fibertype) ;
      

     
     
         if i==1
           pmat1=[E]
         else
            pmat1=[pmat1;E];
         end
             
end

pmat=pmat1';

h = sum(layerthick);

dab = zeros(6,6);


hh = -h/2;

for i=1:n_l
    e1 = pmat(1,i);
    e2 = pmat(2,i);
    e3 = pmat(3,i);

    g12 = pmat(4,i);
    g23 = pmat(5,i);
    g13 = pmat(6,i);

    nu12 = pmat(7,i);
    nu23 = pmat(8,i);
    nu13 = pmat(9,i);
    nu21 = nu12*e2/e1;
    
    c11=e1/(1-nu12*nu21);
    c12=nu12*e2/(1-nu12*nu21);
    c22=e2/(1-nu12*nu21);
    c66 = g12;

    trad = ort(i); % Layer orientation
    n = sind(trad);
    m = cosd(trad);

    m2 = m*m;
    n2 = n*n;
    m2n2 = m2*n2;
    m4 = m2*m2;
    n4 = n2*n2;
    mn3 = m*n2*n;
    m3n = m2*m*n;
    
    %	TRANSFORMATION  (MECHANICAL PART)
    q(1,1) = c11*m4+2*(c12+2*c66)*m2n2+c22*n4;
    q(1,2) = (c11+c22-4*c66)*m2n2+c12*(m4+n4);
    q(1,3) = (c11-c12-2*c66)*m3n+(c12-c22+2*c66)*mn3;
    q(2,2) = c22*m4+2*(c12+2*c66)*m2n2+c11*n4;
    q(2,3) = (c12-c22+2*c66)*m3n+(c11-c12-2*c66)*mn3;
    q(3,3) = (c11-2*c12+c22-2*c66)*m2n2+c66*(m4+n4);
    
    h1 = hh;
    h2 = h1+layerthick(i);
    h12 = h2-h1;
    hh=hh+layerthick(i);
    h122 = (h2*h2-h1*h1)/2.0;
    h123 = (h2^3-h1^3)/3.0;
    
    %   MEMBRANE, COUPLING, BENDING  (MECHANICAL)
    dab(1,1) = dab(1,1)+q(1,1)*h12;
    dab(1,2) = dab(1,2)+q(1,2)*h12;
    dab(1,3) = dab(1,3)+q(1,3)*h12;
    dab(1,4) = dab(1,4)+q(1,1)*h122;
    dab(1,5) = dab(1,5)+q(1,2)*h122;
    dab(1,6) = dab(1,6)+q(1,3)*h122;
    dab(2,2) = dab(2,2)+q(2,2)*h12;
    dab(2,3) = dab(2,3)+q(2,3)*h12;
    dab(2,4) = dab(2,4)+q(1,2)*h122;
    dab(2,5) = dab(2,5)+q(2,2)*h122;
    dab(2,6) = dab(2,6)+q(2,3)*h122;
    dab(3,3) = dab(3,3)+q(3,3)*h12;
    dab(3,4) = dab(3,4)+q(1,3)*h122;
    dab(3,5) = dab(3,5)+q(2,3)*h122;
    dab(3,6) = dab(3,6)+q(3,3)*h122;
    dab(4,4) = dab(4,4)+q(1,1)*h123;
    dab(4,5) = dab(4,5)+q(1,2)*h123;
    dab(4,6) = dab(4,6)+q(1,3)*h123;
    dab(5,5) = dab(5,5)+q(2,2)*h123;
    dab(5,6) = dab(5,6)+q(2,3)*h123;
    dab(6,6) = dab(6,6)+q(3,3)*h123;
end
    
A=[dab(1,1) dab(1,2) dab(1,3)
   dab(1,2) dab(2,2) dab(2,3)
   dab(1,3) dab(2,3) dab(3,3)];
B=inv(A);
E11(w)=1/(B(1,1)*h);
E22(w)=1/(B(2,2)*h);
G12(w)=1/(B(3,3)*h);
V12(w)= -B(1,2)/B(1,1);
V21(w)= -B(1,2)/B(2,2);
       
end
com1=[a',E11',E22',G12',V12',V21']




% reference

format short g
a(1)=0;
b(1)=0;
for w=1:3
    
ort=[ a(w) 0  b(w)]
if (w<3)
a(w+1)=a(w)+45;
b(w+1)=b(w)+45;
else
    
a=a;
b=b;
end
n_l=1;
order=[10];
layerthick=[.005];


for(i=1:n_l)
    
    fibertype=order(i);
    
E_matrix_a  = 3;
E_matrix_t  = 3;
rho_matrix= 1200;
nu_matrix_a = 0.3;
nu_matrix_t = 0.3;
G_matrix=E_matrix_a/(2*(1+nu_matrix_a));
% G_matrix=1.2777;

%Properties of Long Fiber or particulate
E_fiber_1_a   = 230;
E_fiber_1_t   = 8;
rho_fiber_1 = 1760;
nu_fiber_1_a  = 0.256;
nu_fiber_1_t  = 0.256;
G_fiber_1=27.3;
% G_fiber_1=E_fiber__a/(2*(1+nu_fiber_1_a));
length_fiber_1 = 40;
Dia_fiber_1    = 1;    
E_fiber_a   = E_fiber_1_a;
E_fiber_t   = E_fiber_1_t;
nu_fiber_a  = nu_fiber_1_a;
nu_fiber_t  = nu_fiber_1_t;
G_fiber=G_fiber_1;
length_fiber = length_fiber_1;
Dia_fiber    = Dia_fiber_1;
    
[E11_r,E22_r,E33_r,G12_r,G13_r,G23_r,nu_v12_r,nu_v23_r,nu_v31_r,v1_r]=Generic(fibertype,E_matrix_a,E_matrix_t,nu_matrix_a,nu_matrix_t,G_matrix,E_fiber_a,E_fiber_t,nu_fiber_a,nu_fiber_t,G_fiber,length_fiber,Dia_fiber); 
 E= [E11_r*10^9,E22_r*10^9,E33_r*10^9,G12_r*10^9,G13_r*10^9,G23_r*10^9,nu_v12_r,nu_v23_r,nu_v31_r];
        if i==1
           pmat1=E;
        else
           pmat1=[pmat1;E];
        end
             
end

pmat=pmat1';

h = sum(layerthick);

dab = zeros(6,6);


hh = -h/2;

for i=1:n_l
    e1 = pmat(1,i);
    e2 = pmat(2,i);
    e3 = pmat(3,i);

    g12 = pmat(4,i);
    g23 = pmat(5,i);
    g13 = pmat(6,i);

    nu12 = pmat(7,i);
    nu23 = pmat(8,i);
    nu13 = pmat(9,i);
    nu21 = nu12*e2/e1;
    
    c11=e1/(1-nu12*nu21);
    c12=nu12*e2/(1-nu12*nu21);
    c22=e2/(1-nu12*nu21);
    c66 = g12;

    trad = ort(i); % Layer orientation
    n = sind(trad);
    m = cosd(trad);

    m2 = m*m;
    n2 = n*n;
    m2n2 = m2*n2;
    m4 = m2*m2;
    n4 = n2*n2;
    mn3 = m*n2*n;
    m3n = m2*m*n;
    
    %	TRANSFORMATION  (MECHANICAL PART)
    q(1,1) = c11*m4+2*(c12+2*c66)*m2n2+c22*n4;
    q(1,2) = (c11+c22-4*c66)*m2n2+c12*(m4+n4);
    q(1,3) = (c11-c12-2*c66)*m3n+(c12-c22+2*c66)*mn3;
    q(2,2) = c22*m4+2*(c12+2*c66)*m2n2+c11*n4;
    q(2,3) = (c12-c22+2*c66)*m3n+(c11-c12-2*c66)*mn3;
    q(3,3) = (c11-2*c12+c22-2*c66)*m2n2+c66*(m4+n4);
    
    h1 = hh;
    h2 = h1+layerthick(i);
    h12 = h2-h1;
    hh=hh+layerthick(i);
    h122 = (h2*h2-h1*h1)/2.0;
    h123 = (h2^3-h1^3)/3.0;
    
    %   MEMBRANE, COUPLING, BENDING  (MECHANICAL)
    dab(1,1) = dab(1,1)+q(1,1)*h12;
    dab(1,2) = dab(1,2)+q(1,2)*h12;
    dab(1,3) = dab(1,3)+q(1,3)*h12;
    dab(1,4) = dab(1,4)+q(1,1)*h122;
    dab(1,5) = dab(1,5)+q(1,2)*h122;
    dab(1,6) = dab(1,6)+q(1,3)*h122;
    dab(2,2) = dab(2,2)+q(2,2)*h12;
    dab(2,3) = dab(2,3)+q(2,3)*h12;
    dab(2,4) = dab(2,4)+q(1,2)*h122;
    dab(2,5) = dab(2,5)+q(2,2)*h122;
    dab(2,6) = dab(2,6)+q(2,3)*h122;
    dab(3,3) = dab(3,3)+q(3,3)*h12;
    dab(3,4) = dab(3,4)+q(1,3)*h122;
    dab(3,5) = dab(3,5)+q(2,3)*h122;
    dab(3,6) = dab(3,6)+q(3,3)*h122;
    dab(4,4) = dab(4,4)+q(1,1)*h123;
    dab(4,5) = dab(4,5)+q(1,2)*h123;
    dab(4,6) = dab(4,6)+q(1,3)*h123;
    dab(5,5) = dab(5,5)+q(2,2)*h123;
    dab(5,6) = dab(5,6)+q(2,3)*h123;
    dab(6,6) = dab(6,6)+q(3,3)*h123;
end
    
A=[dab(1,1) dab(1,2) dab(1,3)
   dab(1,2) dab(2,2) dab(2,3)
   dab(1,3) dab(2,3) dab(3,3)];
B=inv(A);
E11_R(w)=1/(B(1,1)*h);
E22_R(w)=1/(B(2,2)*h);
G12_R(w)=1/(B(3,3)*h);
V12_R(w)= -B(1,2)/B(1,1);
V21_R(w)= -B(1,2)/B(2,2);
       
end
com2=[a',E11_R',E22_R',G12_R',V12_R',V21_R']


xlswrite('/Users/Abhijeet/Documents/studies/interlaminar shear/Simulations/Laminate.csv',com1,'sheet1'' , A1:F5');
xlswrite('/Users/Abhijeet/Documents/studies/interlaminar shear/Simulations/LaminateRef.csv',com2,'sheet2' , 'H1:M5');

            subplot(3,3,1)
            plot(a,E11,'-.k*',a,E11_R);
            legend('with nano','without nano')
            xlabel('orientation')
            ylabel('Elastic modulus (E axial)');
            
            
            subplot(3,3,3)
            plot(a,E22,'-.k*',a,E22_R);
            legend('with nano','without nano')
            xlabel('orientation')
            ylabel('Elasttic modulus (E transverse_2)');
            
   
            subplot(3,3,5)
            plot(a,G12,'-.rs',a,G12_R);
            legend('with nano','without nano')
            xlabel('orientation')
            ylabel('shear modulus (G12)');
           

            subplot(3,3,7)
            plot(a,V12,'-.cy+',a,V12_R);
            legend('with nano','without nano')
            xlabel('volume fraction (v1)')
            ylabel('Poissons Ratio_v12');
            
            
            subplot(3,3,9)
            plot(a,V21,'-.cy+',a,V21_R);
            legend('with nano','without nano')
            xlabel('volume fraction (v1)')
            ylabel('Poissons Ratio_v21');
          

