%Multi Step homogenization of Nanocomposites using Method of Cells
%Original model developed by Jacob Aboudi
%Developed by Abhinav Alva, Scientist Fellow, STTD, NAL
%Under the Guidance of Dr. S. Raja, Group Head, DAS, STTD, NAL


% STEP1: reference Homogenization for long fiber
% agglomeration=0, computes E_composite for isolated random fibers 
function [E]= MOC_3D_Hybrid(fibertype)

%Properties of Matrix (Isotropic)
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

%Properties of Nanofibers
E_fiber_2_a   = 1000;
E_fiber_2_t   = 1000;
rho_fiber_2 = 2100;
nu_fiber_2_a  = 0.33;
nu_fiber_2_t  = 0.33;
G_fiber_2=E_fiber_2_a/(2*(1+nu_fiber_2_a));
%Fiber Dimensions
length_fiber_2 = 30e-6;
Dia_fiber_2    = 2e-9;
aspect_ratio   = length_fiber_2/Dia_fiber_2;

% %Properties of Particulates
% E_fiber_2_a   = 70;
% E_fiber_2_t   = 70;
% rho_fiber_2 = 2700;
% nu_fiber_2_a  = 0.23;
% nu_fiber_2_t  = 0.23;
% G_fiber_2=E_fiber_2_a/(2*(1+nu_fiber_2_a));
% %Fiber Dimensions
% length_fiber_2 = 1;
% Dia_fiber_2    = 1;
% aspect_ratio   = length_fiber_2/Dia_fiber_2;

%Reinforcement weight Percentage
wt=(2.0/100);




E_fiber_a   = E_fiber_1_a;
E_fiber_t   = E_fiber_1_t;
nu_fiber_a  = nu_fiber_1_a;
nu_fiber_t  = nu_fiber_1_t;
G_fiber=G_fiber_1;

length_fiber = length_fiber_1;
Dia_fiber    = Dia_fiber_1;


% STEP1: First Homogenization
% agglomeration=0, computes E_composite for isolated random fibers 
if (fibertype==7||fibertype==9||fibertype==11)
    agglomeration=3;
else
    agglomeration=0;
end

E_fiber_a   = E_fiber_2_a;
E_fiber_t   = E_fiber_2_t;
rho_fiber = rho_fiber_2;
nu_fiber_a  = nu_fiber_2_a;nu_fiber_t  = nu_fiber_2_t;
G_fiber=G_fiber_2;
length_fiber = length_fiber_2;
Dia_fiber    = Dia_fiber_2;
[E11,E22,E33,G12,G13,G23,nu_v12,nu_v23,nu_v31,l_cell,d_cell,v1]=Nano(fibertype,agglomeration,E_matrix_a,E_matrix_t,rho_matrix,nu_matrix_a,nu_matrix_t,G_matrix,E_fiber_a,E_fiber_t,rho_fiber,nu_fiber_a,nu_fiber_t,G_fiber,wt,length_fiber,Dia_fiber);
E_uagg_a=E11
E_uagg_t=E22;
nu_uagg_a=nu_v12;nu_uagg_t=nu_v23;
if ((fibertype==7)||(fibertype==9)||(fibertype==11))
    G_uagg=G12;
else
    G_uagg=(E_uagg_a/(2*(1+nu_uagg_a)));
end

% STEP2: Second Homogenization
% agglomeration=1, computes E_composite for isolated agglomerated random fibers 
if (fibertype==7||fibertype==9||fibertype==11)
    agglomeration=3;
else
    agglomeration=1;
end
[E11,E22,E33,G12,G13,G23,nu_v12,nu_v23,nu_v31,l_cell,d_cell,v1]=Nano(fibertype,agglomeration,E_matrix_a,E_matrix_t,rho_matrix,nu_matrix_a,nu_matrix_t,G_matrix,E_fiber_a,E_fiber_t,rho_fiber,nu_fiber_a,nu_fiber_t,G_fiber,wt,length_fiber,Dia_fiber);
E_agg_a=E11
E_agg_t=E22;
nu_agg_a=nu_v12;nu_agg_t=nu_v23;
if ((fibertype==7)||(fibertype==9)||(fibertype==11))
    G_agg=G12;
else
    G_agg=(E_agg_a/(2*(1+nu_agg_a)));
end
l_new=l_cell;
d_new=d_cell;

% STEP3: Intermediate Homogenization
% agglomeration=2, computes E_composite for combination of 2 cases stated above
if (fibertype==7||fibertype==9||fibertype==11)
    agglomeration=3;
else
    agglomeration=2;
end
E_matrix_a=E_uagg_a;
E_matrix_t=E_uagg_t;
nu_matrix_a=nu_uagg_a;
nu_matrix_t=nu_uagg_t;
G_matrix=G_uagg;

E_fiber_a=E_agg_a;
E_fiber_t=E_agg_t;
nu_fiber_a=nu_uagg_a;
nu_fiber_t=nu_uagg_t;
G_fiber=G_agg;
length_fiber=d_new;
Dia_fiber=l_new;


[E11,E22,E33,G12,G13,G23,nu_v12,nu_v23,nu_v31,l_cell,d_cell,v1]=Nano(fibertype,agglomeration,E_matrix_a,E_matrix_t,rho_matrix,nu_matrix_a,nu_matrix_t,G_matrix,E_fiber_a,E_fiber_t,rho_fiber,nu_fiber_a,nu_fiber_t,G_fiber,wt,length_fiber,Dia_fiber);


% STEP4: Hybrid Homogenization
% agglomeration=3, computes E_composite for combination of 2 fibers (long & short)
if (fibertype==11)
E_matrix_t=E11;
E_matrix_a=E22;
nu_matrix_t=nu_v12;
nu_matrix_a=nu_v23;
G_matrix=G12;    
else   
E_matrix_a=E11;
E_matrix_t=E22;
nu_matrix_a=nu_v12;
nu_matrix_t=nu_v23;
if ((fibertype==7)||(fibertype==9))
    G_matrix=G12;
else
    G_matrix=(E_matrix_a/(2*(1+nu_matrix_a)));
end
end
E_fiber_a=E_fiber_1_a;
E_fiber_t=E_fiber_1_t;
rho_fiber=rho_fiber_1;
nu_fiber_a=nu_fiber_1_a;
nu_fiber_t=nu_fiber_1_t;
length_fiber=length_fiber_1;
Dia_fiber=Dia_fiber_1;
G_fiber=G_fiber_1;
[E11,E22,E33,G12,G13,G23,nu_v12,nu_v23,nu_v31,v1]=Generic(fibertype,E_matrix_a,E_matrix_t,nu_matrix_a,nu_matrix_t,G_matrix,E_fiber_a,E_fiber_t,nu_fiber_a,nu_fiber_t,G_fiber,length_fiber,Dia_fiber);
E=[E11*10^9,E22*10^9,E33*10^9,G12*10^9,G13*10^9,G23*10^9,nu_v12,nu_v23,nu_v31];

end

            











    
