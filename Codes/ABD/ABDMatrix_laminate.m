


% layerwise E,v,G values
pmat1= [  117e+09   7.28e+09  7.28e+09     3.01e+09     0     0     0.2758    0    0 
          117e+09   7.28e+09  7.28e+09     3.01e+09     0     0     0.2758    0    0
          117e+09   7.28e+09  7.28e+09     3.01e+09     0     0     0.2758    0    0 ] ; 

pmat=pmat1';

           
% layerwise orientation        
ort=[ 90 0  90];

%layerwise thickness in mm
layerthick = [.005 .005 .005];
h = sum(layerthick);

dab = zeros(6,6);

format short g
hh = -h/2;
for i=1:3 %number of iterations for each layer
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
%     nu31 = nu13*e3/e1;
%     nu32 = nu23*e3/e2;
% 
%     rho = pmat(10,i);
    
%     f = 1-nu12*nu21-nu13*nu31-nu23*nu32-nu12*nu23*nu31-nu21*nu13*nu32;

    %	ELASTIC CONSTANTS
%     c11_1 = e1*(1-nu23*nu32)/f;
%     c12_1 = e2*(nu12+nu13*nu32)/f;
%     c13_1 = e3*(nu13+nu12*nu23)/f;
%     c22_1 = e2*(1-nu13*nu31)/f;
%     c23_1 = e3*(nu23+nu21*nu13)/f;
%     c33_1 = e3*(1-nu12*nu21)/f;

    %	ELEMINATION OF eps33  (PLANE STRESS 3D CASE)
    %	MECHANICAL PART (IN-PLANE)
%     c11 = c11_1-(c13_1*(c13_1/c33_1));
%     c12 = c12_1-(c13_1*(c23_1/c33_1));
%     c22 = c22_1-(c23_1*(c23_1/c33_1));
c11=e1/(1-nu12*nu21);
c12=nu12*e2/(1-nu12*nu21);
c22=e2/(1-nu12*nu21);

    %	MECHANICAL PART (SHEAR)
%     c44 = 5/6*g23;
%     c55 = 5/6*g13;
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
    
%     q(4,4) = c44*m2+c55*n2;
%     q(4,5) = (-c44+c55)*m*n;
%     q(5,5) = c44*n2+c55*m2;
%     
    %	POSITION PROPERTIES
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
    
    
%     %   SHEAR  TERMS (MECHANICAL)
%     dab(7,7) = dab(7,7)+q(4,4)*h12;
%     dab(7,8) = dab(7,8)+q(4,5)*h12;
%     dab(8,8) = dab(8,8)+q(5,5)*h12;
end
A=[dab(1,1) dab(1,2) dab(1,3)
   dab(1,2) dab(2,2) dab(2,3)
   dab(1,3) dab(2,3) dab(3,3)]
B=inv(A)

% E,v,G for symmetric cases
E11=1/(B(1,1)*h);
E22=1/(B(2,2)*h);
G12=1/(B(3,3)*h);
V12= -B(1,2)/B(1,1);
V21= -B(1,2)/B(2,2);
C=[E11,E22,G12,V12,V21]
compare=pmat';



