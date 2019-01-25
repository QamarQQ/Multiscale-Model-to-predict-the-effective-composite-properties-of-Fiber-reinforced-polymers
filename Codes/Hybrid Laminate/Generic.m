%Generic Code for Nanocomposites based on Method of Cells by Jacob Aboudi

%Covers 3 Nanoeffects: Equivalent Volume of Nanoparticle, Random
%orientation of fiber & Nanoparticle agglomeration

%Code Developed by Awant Bhagat, BITS-Pilani Student, 2011
%Code Developed by Abhinav Alva, Scientist Fellow, STTD, NAL, 2011-12
%Under the Guidance of Dr. S. Raja, Group Head, DAS, STTD, NAL


 
function[E11,E22,E33,G12,G13,G23,nu_v12,nu_v23,nu_v31,v1]=Generic(fibertype,E_matrix_a,E_matrix_t,nu_matrix_a,nu_matrix_t,G_matrix,E_fiber_a,E_fiber_t,nu_fiber_a,nu_fiber_t,G_fiber,length_fiber,Dia_fiber);

if fibertype==10 
% Basic Characteristic Values for the Matrix are assigned from main code
E_a_m = E_matrix_a;
E_t_m = E_matrix_t;
v_t_m = nu_matrix_t;
G_t_m = G_matrix;
G_a_m = G_matrix;
v_a_m = nu_matrix_a;
k_m = (0.25 * E_a_m)/(0.5*(1 - v_t_m)*(E_a_m/E_t_m) - (v_a_m^2));
else
%Properties of Matrix (Isotropic)
E_matrix_a  = 3;
E_matrix_t  = 3;
rho_matrix= 1200;
nu_matrix_a = 0.3;
nu_matrix_t = 0.3;
G_matrix=E_matrix_a/(2*(1+nu_matrix_a));
% G_matrix=1.2777;
end
% Corresponding Stiffness Matrix Calculation for Matrix
c11_m = E_a_m + 4*k_m*(v_a_m^2);
c12_m = 2*k_m*(v_a_m);
c22_m = k_m + G_t_m;
c23_m = k_m - G_t_m;
c66_m = G_a_m;

% Basic Characteristic Values for the Fibre are assigned from main code
E_a_f = E_fiber_a;
E_t_f = E_fiber_t;
v_t_f = nu_fiber_t;
v_a_f = nu_fiber_a;
G_t_f = G_fiber;
G_a_f = G_fiber;
k_f = (0.25 * E_a_f)/(0.5*(1 - v_t_f)*(E_a_f/E_t_f) - (v_a_f^2));

% Corresponding Stiffness Matrix Calculation for Fibre
c11_f = E_a_f + 4*k_f*(v_a_f^2);
c12_f = 2*k_f*(v_a_f);
c22_f = k_f + G_t_f;
c23_f = k_f - G_t_f;
c66_f = G_a_f;

v_m=1;
vol_fraction = .502;

            if (fibertype==1||fibertype==9||fibertype==10||fibertype==11)
            d1=1e15;
            l1=1;
            h1=l1;
            a_r = (d1/l1);
            elseif (fibertype==2)
            d1=1e15;
            l1=1e15;
            h1=1;
            a_r = (d1/l1);
            elseif (fibertype==3||fibertype==7||fibertype==8)
            d1=1;
            l1=1;
            h1=1;
            a_r = (d1/l1);
            else
            d1=length_fiber;
            l1=Dia_fiber;
            h1=l1;
            a_r = (d1/l1);
            end
            
    for vol_it=1:v_m
                    
            syms xyz a1;
            a1 = solve(vol_fraction * (xyz*xyz+2*xyz+1)* (a_r*xyz + 1) -(a_r*xyz*xyz*xyz), xyz);
            a1 = real (eval(a1(1)));
            h2 = h1/a1;
            l2 = h2;
            d2 = h2;
            d = (d1 + d2);
            h = h1 + h2;
            l = l1 + l2;
                       
            %Calculation of G12 Value
            
            num1 = (d1*h1*c66_m + d2*h2*c66_f + d2*h1*c66_f + d1*h2*c66_f);
            N111 = d*h*c66_m/num1;
            N121 = d*h*c66_f/num1;
            N211 = d*h*c66_f/num1;
            N221 = d*h*c66_f/num1;

            num2 = d1*h1 + d2*h2 + d2*h1 + d1*h2;
            N122 = d*h/num2;
            N212 = d*h/num2;
            N222 = d*h/num2;
            N112 = d*h/num2;

            G12_1 = (1/(d*h*l))  *  (d1*h1*l1*c66_f*N111 + d1*h1*l2*N112*c66_m + d1*h2*l1*c66_m*N121 + d1*h2*l2*c66_m*N122 + d2*h1*l1*c66_m*N211 + d2*h1*l2*c66_m*N212 + d2*h2*l1*c66_m*N221 + d2*h2*l2*c66_m*N222);
            
                      
            num3 = (d1*l1*c66_m + d2*l2*c66_f + d2*l1*c66_f + d1*l2*c66_f);
            N111 = l*d*c66_m/num3;
            N112 = l*d*c66_f/num3;
            N211 = l*d*c66_f/num3;
            N212 = l*d*c66_f/num3;

            num4 = l1*d1 + l2*d2 + l2*d1 + l1*d2;
            N122 = l*d/num4;
            N121 = l*d/num4;
            N221 = l*d/num4;
            N222 = l*d/num4;

            G13_1 = (1/(d*h*l))  *  (d1*h1*l1*c66_f*N111 + d1*h1*l2*N112*c66_m + d1*h2*l1*c66_m*N121 + d1*h2*l2*c66_m*N122 + d2*h1*l1*c66_m*N211 + d2*h1*l2*c66_m*N212 + d2*h2*l1*c66_m*N221 + d2*h2*l2*c66_m*N222);

            %Calculation of G23 Value
            
            c_eff_fiber  = ((c22_f - c23_f)/2);
            c_eff_matrix = ((c22_m - c23_m)/2);

            num5 = (l1*h1*c_eff_matrix + l2*h1*c_eff_fiber + l1*h2*c_eff_fiber + l2*h2*c_eff_fiber);
            N111 = l*h*c_eff_matrix/num5;
            N122 = l*h*c_eff_fiber/num5;
            N112 = l*h*c_eff_fiber/num5;
            N121 = l*h*c_eff_fiber/num5;

            num6 = l1*h1 + l2*h2 + l2*h1 + l1*h2;
            N211 = l*h/num6;
            N212 = l*h/num6;
            N222 = l*h/num6;
            N221 = l*h/num6;

            G23_1 = (1/(d*h*l))  *  (d1*h1*l1*c_eff_fiber*N111 + d1*h1*l2*N112*c_eff_matrix + d1*h2*l1*c_eff_matrix*N121 + d1*h2*l2*c_eff_matrix*N122 + d2*h1*l1*c_eff_matrix*N211 + d2*h1*l2*c_eff_matrix*N212 + d2*h2*l1*c_eff_matrix*N221 + d2*h2*l2*c_eff_matrix*N222);
                        
            % G Calculation Over
            % E Calculation Started
            % Relation between elementary values and A's
            
            A1  = c11_f + c11_m*(d1/d2); 
            A2  = c12_f;
            A3  = c12_f;
            A4  = c12_m *(h2/h1);
            A5  = c12_m;
            A6  = c11_m + c11_m*(d1/d2); 
            A7  = c12_m;
            A8  = c12_m *(h2/h1);
            A9  = c12_m*(l1/l2);
            A10 = c12_m*(l1/l2);
            A11 = c11_m + c11_m*(d2/d1);
            A12 = c12_m*(h1/h2);
            A13 = c12_m*(l2/l1);
            A14 = c12_m;
            A15 = c12_m*(l2/l1);
            A16 = c11_m + c11_m *(d2/d1);
            A17 = c12_m;
            A18 = c12_m*(h1/h2);
            A19 = c12_m;
            A20 = c12_m;
            A21 = c12_f;
            A22 = c22_f + c22_m*(h1/h2);
            A23 = c23_f;
            A24 = c12_m*(d2/d1);
            A25 = c23_m*(l2/l1);
            A26 = c12_m;
            A27 = c22_m + c22_m*(h1/h2);
            A28 = c12_m*(d2/d1);
            A29 = c23_m*(l1/l2);
            A30 = c23_m;%VERIFIED
            A31 = c12_m;
            A32 = c22_m + c22_m*(h2/h1);
            A33 = c12_m*(d1/d2);
            A34 = c23_m*(l2/l1);
            A35 = c23_m;
            A36 = c12_m;
            A37 = c22_m + c22_m*(h2/h1);
            A38 = c23_m;
            A39 = c23_m*(l1/l2);
            A40 = c12_m*(d1/d2);
            A41 = c12_f;
            A42 = c23_f;
            A43 = c22_f + c22_m*(l1/l2);
            A44 = c12_m;
            A45 = c23_m;
            A46 = c12_m*(d2/d1);
            A47 = c23_m*(h1/h2);
            A48 = c22_m + c22_m*(l2/l1);
            A49 = c23_m*(h1/h2);
            A50 = c12_m*(d2/d1);
            A51 = c22_m + c22_m*(l1/l2);
            A52 = c12_m*(d1/d2);
            A53 = c23_m*(h2/h1);
            A54 = c12_m*(d1/d2);
            A55 = c23_m*(h2/h1);
            A56 = c12_m;
            A57 = c23_m;
            A58 = c22_m + c22_m*(l2/l1);
            A59 = c12_m;
            A60 = c23_m;
            
            % The Matrix consisting of the above to find the values of phi, sy and chi.
            
            A   = [A1 0 0 0 A2 0 A4 0 A3 0 -A5 0;
                0 A6 0 0 0 A7 0 A8 -A10 0 A9 0;
                0 0 A11 0 A12 0 A14 0 0 A13 0 -A15;
                0 0 0 A16 0 A18 0 A17 0 -A20 0 A19;
                A21 0 A24 0 A22 0 0 0 A23 A25 0 0;
                0 A26 0 A28 0 A27 0 0 -A29 -A30 0 0;
                A33 0 A31 0 0 0 A32 0 0 0 -A35 -A34;
                0 A40 0 A36 0 0 0 A37 0 0 A39 A38;
                A41 -A44 0 0 A42 -A45 0 0 A43 0 0 0;
                0 0 A46 -A50 A47 -A49 0 0 0 A48 0 0;
                -A54 A52 0 0 0 0 -A55 A53 0 0 A51 0;
                0 0 -A59 A56 0 0 -A60 A57 0 0 0 A58];

            syms e11 e22 e33
            
            % Relation between elementary values and J's
            
            J1  = vpa (c11_m *(d/d2)*e11 + c12_m*(h/h1)*e22);
            J2  = vpa (c11_m *(d/d2)*e11 + c12_m*(h/h1)*e22);
            J3  = vpa (c11_m *(d/d1)*e11 + c12_m*(h/h2)*e22); 
            J4  = vpa (c11_m *(d/d1)*e11 + c12_m*(h/h2)*e22);
            J5  = vpa (c12_m *(d/d1)*e11 + c22_m*(h/h2)*e22 + c23_m*(l/l1)*e33);
            J6  = vpa (c12_m *(d/d1)*e11 + c22_m*(h/h2)*e22 - c23_m*(l/l2)*e33);
            J7  = vpa (c12_m *(d/d2)*e11 + c22_m*(h/h1)*e22 - c23_m*(l/l1)*e33);
            J8  = vpa (c12_m *(d/d2)*e11 + c22_m*(h/h1)*e22 + c23_m*(l/l2)*e33);
            J9  = vpa (c22_m*(l/l2)*e33);
            J10 = vpa (c22_m*(l/l1)*e33);
            J11 = vpa (c22_m*(l/l2)*e33);
            J12 = vpa (c22_m*(l/l1)*e33);

            J=[ J1; J2; J3; J4; J5; J6; J7; J8; J9; J10; J11; J12];

            Result  = A\J; %Solving through Gaussian Elimination

            phi1111 = Result (1,1);
            phi1112 = Result (2,1);
            phi1221 = Result (3,1);
            phi1222 = Result (4,1);
            sy2111  = Result (5,1);
            sy2112  = Result (6,1);
            sy2221  = Result (7,1);
            sy2222  = Result (8,1);
            chi3111 = Result (9,1);
            chi3122 = Result (10,1);
            chi3211 = Result (11,1);
            chi3222 = Result (12,1);

            phi1211 = (d*e11 - d1*phi1111)/d2;
            phi1212 = (d*e11 - d1*phi1112)/d2;
            phi1121 = (d*e11 - d2*phi1221)/d1;
            phi1122 = (d*e11 - d2*phi1222)/d1;
            sy2121  = (h*e22 - h1*sy2111)/h2;
            sy2122  = (h*e22 - h1*sy2112)/h2;
            sy2211  = (h*e22 - h2*sy2221)/h1;
            sy2212  = (h*e22 - h2*sy2222)/h1;
            chi3112 = (l*e33 - l1*chi3111)/l2;
            chi3121 = (l*e33 - l2*chi3122)/l1;
            chi3212 = (l*e33 - l1*chi3211)/l2;
            chi3221 = (l*e33 - l2*chi3222)/l1;

            S11111  = c11_f*phi1111 + c12_f*(sy2111 + chi3111);  
            S11112  = c11_m*phi1112 + c12_m*(sy2112 + chi3112);
            S11121  = c11_m*phi1121 + c12_m*(sy2121 + chi3121);
            S11122  = c11_m*phi1122 + c12_m*(sy2122 + chi3122);
            S11211  = c11_m*phi1211 + c12_m*(sy2211 + chi3211);
            S11221  = c11_m*phi1221 + c12_m*(sy2221 + chi3221);
            S11212  = c11_m*phi1212 + c12_m*(sy2212 + chi3212);
            S11222  = c11_m*phi1222 + c12_m*(sy2222 + chi3222);

            S22111  = c12_f*phi1111 + c22_f*sy2111 + c23_f*chi3111;
            S22112  = c12_m*phi1112 + c22_m*sy2112 + c23_m*chi3112;
            S22121  = c12_m*phi1121 + c22_m*sy2121 + c23_m*chi3121;
            S22122  = c12_m*phi1122 + c22_m*sy2122 + c23_m*chi3122;
            S22211  = c12_m*phi1211 + c22_m*sy2211 + c23_m*chi3211;
            S22221  = c12_m*phi1221 + c22_m*sy2221 + c23_m*chi3221;
            S22212  = c12_m*phi1212 + c22_m*sy2212 + c23_m*chi3212;
            S22222  = c12_m*phi1222 + c22_m*sy2222 + c23_m*chi3222;

            S33111  = c12_f*phi1111 + c23_f*sy2111 + c22_f*chi3111;
            S33112  = c12_m*phi1112 + c23_m*sy2112 + c22_m*chi3112;
            S33121  = c12_m*phi1121 + c23_m*sy2121 + c22_m*chi3121;
            S33122  = c12_m*phi1122 + c23_m*sy2122 + c22_m*chi3122;
            S33211  = c12_m*phi1211 + c23_m*sy2211 + c22_m*chi3211;
            S33221  = c12_m*phi1221 + c23_m*sy2221 + c22_m*chi3221;
            S33212  = c12_m*phi1212 + c23_m*sy2212 + c22_m*chi3212;
            S33222  = c12_m*phi1222 + c23_m*sy2222 + c22_m*chi3222;

            % Normal Stress sig11 and final value of E11
            
            clear a b sol_1 sol_2 sol_11 sol_22 tempe11 tempe22 tempe33 e11 e22 e33
            syms e11 e22 e33

            a = (1/(d*h*l))*(d1*h1*l1*S22111 + d1*h1*l2*S22112 + d1*h2*l1*S22121 + d1*h2*l2*S22122 + d2*h1*l1*S22211 + d2*h2*l1*S22221 + d2*h1*l2*S22212 + d2*h2*l2*S22222);
            b = (1/(d*h*l))*(d1*h1*l1*S33111 + d1*h1*l2*S33112 + d1*h2*l1*S33121 + d1*h2*l2*S33122 + d2*h1*l1*S33211 + d2*h2*l1*S33221 + d2*h1*l2*S33212 + d2*h2*l2*S33222);
            sol_1 = solve(a,e22);
            sol_2 = solve(b,e22);
            tempe33 = solve(sol_1 - sol_2,e33);

            a = (1/(d*h*l))*(d1*h1*l1*S22111 + d1*h1*l2*S22112 + d1*h2*l1*S22121 + d1*h2*l2*S22122 + d2*h1*l1*S22211 + d2*h2*l1*S22221 + d2*h1*l2*S22212 + d2*h2*l2*S22222);
            b = (1/(d*h*l))*(d1*h1*l1*S33111 + d1*h1*l2*S33112 + d1*h2*l1*S33121 + d1*h2*l2*S33122 + d2*h1*l1*S33211 + d2*h2*l1*S33221 + d2*h1*l2*S33212 + d2*h2*l2*S33222);
            sol_11 = solve(a,e33);
            sol_22 = solve(b,e33);
            tempe22 = solve(sol_11 - sol_22,e22);

            e22 = tempe22;
            e33 = tempe33;

            sig11   = (1/(d*h*l))*(d1*h1*l1*S11111 + d1*h1*l2*S11112 + d1*h2*l1*S11121 + d1*h2*l2*S11122 + d2*h1*l1*S11211 + d2*h2*l1*S11221 + d2*h1*l2*S11212 + d2*h2*l2*S11222);

            E11_1  = eval(sig11/e11);
            Poisson_Ratio_v12_1 = eval(-1*(e22/e11));
            
            
             %Normal Stress sig22 and final value of E22
            
                clear a b sol_1 sol_2 sol_11 sol_22 tempe11 tempe22 tempe33 e11 e22 e33
                syms e11 e22 e33

                a = (1/(d*h*l))*(d1*h1*l1*S11111 + d1*h1*l2*S11112 + d1*h2*l1*S11121 + d1*h2*l2*S11122 + d2*h1*l1*S11211 + d2*h2*l1*S11221 + d2*h1*l2*S11212 + d2*h2*l2*S11222);
                b = (1/(d*h*l))*(d1*h1*l1*S33111 + d1*h1*l2*S33112 + d1*h2*l1*S33121 + d1*h2*l2*S33122 + d2*h1*l1*S33211 + d2*h2*l1*S33221 + d2*h1*l2*S33212 + d2*h2*l2*S33222);
                sol_1 = solve(a,e33);
                sol_2 = solve(b,e33);
                tempe11 = solve(sol_1 - sol_2,e11);

                a = (1/(d*h*l))*(d1*h1*l1*S11111 + d1*h1*l2*S11112 + d1*h2*l1*S11121 + d1*h2*l2*S11122 + d2*h1*l1*S11211 + d2*h2*l1*S11221 + d2*h1*l2*S11212 + d2*h2*l2*S11222);
                b = (1/(d*h*l))*(d1*h1*l1*S33111 + d1*h1*l2*S33112 + d1*h2*l1*S33121 + d1*h2*l2*S33122 + d2*h1*l1*S33211 + d2*h2*l1*S33221 + d2*h1*l2*S33212 + d2*h2*l2*S33222);
                sol_11 = solve(a,e11);
                sol_22 = solve(b,e11);
                tempe33 = solve(sol_11 - sol_22,e33);

                e11 = tempe11;
                e33 = tempe33;

                sig22 = (1/(d*h*l))*(d1*h1*l1*S22111 + d1*h1*l2*S22112 + d1*h2*l1*S22121 + d1*h2*l2*S22122 + d2*h1*l1*S22211 + d2*h2*l1*S22221 + d2*h1*l2*S22212 + d2*h2*l2*S22222);
                E22_1   = eval(sig22/e22);
                Poisson_Ratio_v23_1 = eval(-1*(e33/e22));

                %Normal Stress sig33 and final value of E33
                clear a b sol_1 sol_2 sol_11 sol_22 tempe11 tempe22 tempe33 e11 e22 e33
                syms e11 e22 e33

                a = (1/(d*h*l))*(d1*h1*l1*S11111 + d1*h1*l2*S11112 + d1*h2*l1*S11121 + d1*h2*l2*S11122 + d2*h1*l1*S11211 + d2*h2*l1*S11221 + d2*h1*l2*S11212 + d2*h2*l2*S11222);
                b = (1/(d*h*l))*(d1*h1*l1*S22111 + d1*h1*l2*S22112 + d1*h2*l1*S22121 + d1*h2*l2*S22122 + d2*h1*l1*S22211 + d2*h2*l1*S22221 + d2*h1*l2*S22212 + d2*h2*l2*S22222);
                sol_1 = solve(a,e22);
                sol_2 = solve(b,e22);
                tempe11 = solve(sol_1 - sol_2,e11);
                a = (1/(d*h*l))*(d1*h1*l1*S11111 + d1*h1*l2*S11112 + d1*h2*l1*S11121 + d1*h2*l2*S11122 + d2*h1*l1*S11211 + d2*h2*l1*S11221 + d2*h1*l2*S11212 + d2*h2*l2*S11222);
                b = (1/(d*h*l))*(d1*h1*l1*S22111 + d1*h1*l2*S22112 + d1*h2*l1*S22121 + d1*h2*l2*S22122 + d2*h1*l1*S22211 + d2*h2*l1*S22221 + d2*h1*l2*S22212 + d2*h2*l2*S22222);
                sol_11 = solve(a,e11);
                sol_22 = solve(b,e11);
                tempe22 = solve(sol_11 - sol_22,e22);

                e22 = tempe22;
                e11 = tempe11;

                sig33   = (1/(d*h*l))*(d1*h1*l1*S33111 + d1*h1*l2*S33112 + d1*h2*l1*S33121 + d1*h2*l2*S33122 + d2*h1*l1*S33211 + d2*h2*l1*S33221 + d2*h1*l2*S33212 + d2*h2*l2*S33222);
                E33_1     = eval(sig33/e33);
                Poisson_Ratio_v31_1 = eval(-1*(e11/e33));
            
           
            
            v1(vol_it)=vol_fraction;
            vol_fraction=vol_fraction+0.1;
            G12(vol_it)=G12_1;
            G13(vol_it)=G13_1;
            G23(vol_it)=G23_1;
            E11(vol_it)=E11_1;
            E22(vol_it)=E22_1;
            E33(vol_it)=E33_1;
            nu_v12(vol_it) = Poisson_Ratio_v12_1;
            nu_v23(vol_it) = Poisson_Ratio_v23_1;
            nu_v31(vol_it) = Poisson_Ratio_v31_1;
                         
                         
    end
    v1;
    G12;G13;G23;
    E11=eval(E11);
    E22=eval(E22);
    E33=eval(E33);
    nu_v12;nu_v23;nu_v31;
    
               




