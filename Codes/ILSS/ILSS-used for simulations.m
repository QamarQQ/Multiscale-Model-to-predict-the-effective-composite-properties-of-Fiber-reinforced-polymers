%INPUTS

N=3;


%Txz calcuation
%E1=117/(1-(.2758*.017141));    % c {0 }
%E1=8.44/(1-(.403531*.403531)); % c {45}
%E2=7.28/(1-(.017141*.2758));   % c {90}

%Tyz calculation
%E1=7.28/(1-(.017141*.2758));    % c {0 }
%E1=8.44/(1-(.403531*.403531));  % c {45}
%E2=117/(1-(.2758*.017141));     % c {90}


% weight percentage=1

  %Txz calcuation 
%E2=117.06/(1-(.31972*.03994)); % t {0 } %E2=119.26/(1-(.2767*.01697)); % a {0 }  %E2=118.24/(1-(.2779*.02625)); % r {0 } 
%E2=9.5981/(1-(.65707*.65707)); % t {45} %E2=8.2353/(1-(.4128*.4128)); % a {45}   %E2=14.32/(1-(.35622*.35622)); % r {45} 
%E2=14.625/(1-(.31972*.03994)); % t {90} %E2=7.32/(1-(.2767*.01697)); % a {90}    %E2=11.172/(1-(.2779*.02625)); % r {90}
  
  %Tyz calculation
%E2=14.625/(1-(.31972*.03994)); % t {0 } %E2=7.32/(1-(.2767*.01697));   % a {0 }   %E2=11.172/(1-(.2779*.02625)); % r {0 } 
%E2=9.5981/(1-(.65707*.65707)); % t {45} %E2=8.2353/(1-(.4128*.4128));  % a {45}   %E2=14.32/(1-(.35622*.35622)); % r {45} 
%E2=117.06/(1-(.31972*.03994)); % t {90} %E2=119.26/(1-(.2767*.01697)); % a {90}   %E2=118.24/(1-(.2779*.02625)); % r {90}  
 
  
  
%weight percentage=2

%Tyz calculation
   %E2=15.35/(1-(.2796*.0357)); % r {0 } 
   %E2=18.897/(1-(.31425*.31425)); % r {45} 
   %E2=120.21/(1-(.2796*.0357)); % r {90}


 
  E1=15.35/(1-(.2796*.0357)); % r {0 }    
  E2=120.21/(1-(.2796*.0357)); % r {90} 
  E3=15.35/(1-(.2796*.0357)); % r {0 }    
%  
  E= [E1 E2 E3 ]
%E=[80  800 87 ]
%E=[10 20 30 40];
%E=[ (25/(1-(.25*.01))) (1/(1-(.25*.01))) (25/(1-(.2r5*.01)))  (1/(1-(.25*.01)))  (25/(1-(.25*.01)))   ];    % 0 90 0 90 0
%E=[ (25/(1-(.25*.01))) (1/(1-(.25*.01))) (1.325/(1-(.3245)^2))  (1/(1-(.25*.01)))  (25/(1-(.25*.01)))   ]  % 0 90 45 90 0 
%E= [ (2.192/(1-(.4082*.2994  ))) (1/(1-(.25*.01))) (1.325/(1-(.3245)^2)) (1.068/(1-(.1989*.40822))) 2.4971]; % -30 90 45 60 30
% E=[40 10]
b= 10;
t=[.5  .5  .5];
Q=1;
%Calculation of boundaries of layers : z_1(i) or zk

c=sum(t)/2;
z_1(1)=c;
for i=1:N
    
    z_1(i+1)=z_1(i)-(t(i))
end

clear i;

%Calculation of normalization factors : n(i) 

E_a=sum(E)/N;
for i=1:N
    n(i)=E(i)/E_a;
end

clear i;

%Calculation of nuetral axes : z_g or z*
for i=1:N
    c(i)=((z_1(i+1)+z_1(i))*n(i)*t(i));
end

clear i;

a=sum(c)/2;
for i=1:N
    C1(i)=(t(i)*n(i));
end
w=sum(C1);
z_g = a/w;


%Calculation of Boundaries wrt Nuetral axes or z or z_cap

 for i=1:N+1;
 z(i) = z_1(i)-z_g;
end

clear i;

%Calculation of Moment of Inertia

for i=1:N
e(i)= (b/12)* n(i)*(t(i)*t(i)*t(i));
f(i)= (b/4)*n(i)*(t(i))*((z(i+1)+z(i))^2);
end

I = sum(e)+sum(f);

clear i;


% Layerwise/lamina wise shear stress distribution Calculation.

for i=1:N
    
    
    t_21=0;
    t_22=0;
    t_23=0;
    t_2_1=0;
    t_2_2=0;
    t_2_3=0;
    
    r=((z_1(i)-(z_1(i+1)))/100);
    x_L(:,:,i) = [ z_1(i) : -r :(z_1(i+1))]
 
    
    if z_1(i)*z_1(i+1)>=0 
      
      % t-1 calutaion starts  
        if z_1(i)<0
            h1=-1*z_1(i);
        else
        h1=z_1(i);
        end
    
        if z_1(i+1)<0
           h2=-1*z_1(i+1);
        else
           h2=z_1(i+1);
        end 
    
        if h1>h2
          H= z(i)
        else 
         H= z(i+1) 
        end
        t_1(:,:,i) = (b/2)*n(i)*((x_L(:,:, i )-z_g).^2-(H)^2)
    
     % t-1 calutaion ends
     
     % t-2 calutaion starts 
     
                                     if z_1(i)*z_1(i+1)>0 && z_1(i)>0
                                     display case1
                                     k=i-1;
                                     m=k+1;
                                     
                                     
                                         
                                     for m=k+1:N
                                     t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     
                                     
                                     t_21=t_2_1;
                                     t_2_1=0;
                                     t_2(i)=t_21
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x_L(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                  
                                  elseif z_1(i)*z_1(i+1)>0 && z_1(i)<0
                                      display case2
                                      
             
                                     k=i;
                                     m=k+1;
                                     
                                     
                                     for m=k+1:N
                                     t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     t_22=t_2_2;
                                     t_2_2=0;
                                     t_2(i)=t_22
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x_L(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                   
                                  elseif  z_1(i)*z_1(i+1)==0 && z_1(i)==0
                                      display case3
                                      k=i;
                                     m=k+1;
                                     
                                     
                                     for m=k+1:N
                                     t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     t_22=t_2_2;
                                     t_2_2=0;
                                     t_2(i)=t_22
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x_L(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                   
                                  elseif z_1(i)*z_1(i+1)==0 && z_1(i+1)==0
                                     display case4
                                     k=i-1;
                                     m=k+1;
                                     
                                     
                                         
                                     for m=k+1:N
                                     t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     
                                     
                                     t_21=t_2_1;
                                     t_2_1=0;
                                     t_2(i)=t_21
                                     Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                     T_x_L(:,:,i) = -(Q_x(:,:,i)*Q)/I
     
                                     end
     
                                     
         % t-2 calutaion ends                            
     
    else
        
        % t-1 and t-2 calutaion starts
        for v= 1:length(x_L)
           
         if x_L(:,v,i)>=0
         display case6  
         k=i-1;
         m=k+1;
         t_1(:,v,i) = (b/2)*n(i)*((x_L(:,v, i )-z_g).^2-(z(m))^2)
                                        
                                     
                                     
                                         
                                        for m=k+1:N
                                        t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                        end
                                     
                                     
                                        t_21=t_2_1
                                        t_2_1=0;
                                        t_2(i)=t_21
                                        Q_x(:,v,i) = t_1(:,v,i) + t_2(i); 
    
                                        T_x_L(:,v,i) = -(Q_x(:,v,i)*Q)/I
         
         else
         display case7
         k=i;
         m=k+1;
         t_1(:,v,i) = (b/2)*n(i)*((x_L(:,v, i )-z_g).^2-(z(m))^2)
                                        k=i;
                                        m=k+1;
                                     
                                     
                                        for m=k+1:N
                                        t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                        end
                                        t_22=t_2_2
                                        t_2_2=0;
                                        t_2(i)=t_22
                                        Q_x(:,v,i) = t_1(:,v,i) + t_2(i); 
    
                                        T_x_L(:,v,i) = -(Q_x(:,v,i)*Q)/I
       
         end
         
        v=v+1; 
        
        end
        
        
       % t-1and t-2  calutaion ends
    end  
     
                                     
                                   
end



M_L=max(T_x_L(:,:,:))

% PLOTTING %   

% figure (1)
% plot (  T_x(:,:,1), x_L(:,:, 1) ) ;
 
% figure (1)
% plot (  T_x(:,:,1), x_L(:,:, 1) , T_x(:,:,2), x_L(:,:, 2) );
 
%  
%  figure (1)
%  plot (  T_x_L(:,:,1), x_L(:,:, 1) , T_x_L(:,:,2), x_L(:,:, 2) ,T_x_L(:,:,3) , x_L(:,:, 3) );
  

% figure (2)
% plot (  T_x(:,:,1), x_L(:,:, 1) , T_x(:,:,2), x_L(:,:, 2) ,T_x(:,:,3) ,x_L(:,:, 3) , T_x(:,:,4) ,x_L(:,:, 4) ) ;

% figure (2)
% plot (  T_x(:,:,1), x_L(:,:, 1) , T_x(:,:,2), x_L(:,:, 2) ,T_x(:,:,3) ,x_L(:,:, 3) , T_x(:,:,4) ,x_L(:,:, 4) ) ;
 
% figure (2)
% plot (  T_x(:,:,1), x_L(:,:, 1) , T_x(:,:,2), x_L(:,:, 2) ,T_x(:,:,3) ,x_L(:,:, 3) , T_x(:,:,4) ,x_L(:,:, 4)  ,T_x(:,:,5) , x_L(:,:, 5) ) ;

% figure (1)
% plot (  T_x(:,:,1), x_L(:,:, 1) , T_x(:,:,2), x_L(:,:, 2) ,T_x(:,:,3) ,x_L(:,:, 3) , T_x(:,:,4) ,x_L(:,:, 4)  ,T_x(:,:,5) , x_L(:,:, 5) ,T_x(:,:,6) , x_L(:,:, 6) ) ;






































%INPUTS

N=3;

%Txz calcuation
%E1=117/(1-(.2758*.017141));    % c {0 }
%E1=8.44/(1-(.403531*.403531)); % c {45}
%E2=7.28/(1-(.017141*.2758));   % c {90}

%Tyz calculation
%E1=7.28/(1-(.017141*.2758));    % c {0 }
%E1=8.44/(1-(.403531*.403531));  % c {45}
%E2=117/(1-(.2758*.017141));     % c {90}




% 
  E1=7.28/(1-(.017141*.2758));    % c {0 }
  E2=117/(1-(.2758*.017141));     % c {90}
  E3=7.28/(1-(.017141*.2758));    % c {0 }
  
% %  
 E= [E1 E2 E3 ]
% E=[8  8 8 ]
%E=[10 20 30 40];
%E=[ (25/(1-(.25*.01))) (1/(1-(.25*.01))) (25/(1-(.2r5*.01)))  (1/(1-(.25*.01)))  (25/(1-(.25*.01)))   ];    % 0 90 0 90 0
%E=[ (25/(1-(.25*.01))) (1/(1-(.25*.01))) (1.325/(1-(.3245)^2))  (1/(1-(.25*.01)))  (25/(1-(.25*.01)))   ]  % 0 90 45 90 0 
%E= [ (2.192/(1-(.4082*.2994  ))) (1/(1-(.25*.01))) (1.325/(1-(.3245)^2)) (1.068/(1-(.1989*.40822))) 2.4971]; % -30 90 45 60 30
% E=[40 10]
b= 10;
t=[.5  .5  .5];
Q=1;
%Calculation of boundaries of layers : z_1(i) or zk

c=sum(t)/2;
z_1(1)=c;
for i=1:N
    
    z_1(i+1)=z_1(i)-(t(i))
end

clear i;

%Calculation of normalization factors : n(i) 

E_a=sum(E)/N;
for i=1:N
    n(i)=E(i)/E_a;
end

clear i;

%Calculation of nuetral axes : z_g or z*
for i=1:N
    c(i)=((z_1(i+1)+z_1(i))*n(i)*t(i));
end

clear i;

a=sum(c)/2;
for i=1:N
    C1(i)=(t(i)*n(i));
end
w=sum(C1);
z_g = a/w;


%Calculation of Boundaries wrt Nuetral axes or z or z_cap

 for i=1:N+1;
 z(i) = z_1(i)-z_g;
end

clear i;

%Calculation of Moment of Inertia

for i=1:N
e(i)= (b/12)* n(i)*(t(i)*t(i)*t(i));
f(i)= (b/4)*n(i)*(t(i))*((z(i+1)+z(i))^2);
end

I = sum(e)+sum(f);

clear i;


% Layerwise/lamina wise shear stress distribution Calculation.

for i=1:N
    
    
    t_21=0;
    t_22=0;
    t_23=0;
    t_2_1=0;
    t_2_2=0;
    t_2_3=0;
    
    r=((z_1(i)-(z_1(i+1)))/100);
    x_R(:,:,i) = [ z_1(i) : -r :(z_1(i+1))]
 
    
    if z_1(i)*z_1(i+1)>=0 
      
      % t-1 calutaion starts  
        if z_1(i)<0
            h1=-1*z_1(i);
        else
        h1=z_1(i);
        end
    
        if z_1(i+1)<0
           h2=-1*z_1(i+1);
        else
           h2=z_1(i+1);
        end 
    
        if h1>h2
          H= z(i)
        else 
         H= z(i+1) 
        end
        t_1(:,:,i) = (b/2)*n(i)*((x_R(:,:, i )-z_g).^2-(H)^2)
    
     % t-1 calutaion ends
     
     % t-2 calutaion starts 
     
                                     if z_1(i)*z_1(i+1)>0 && z_1(i)>0
                                     display case1
                                     k=i-1;
                                     m=k+1;
                                     
                                     
                                         
                                     for m=k+1:N
                                     t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     
                                     
                                     t_21=t_2_1;
                                     t_2_1=0;
                                     t_2(i)=t_21
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x_R(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                  
                                  elseif z_1(i)*z_1(i+1)>0 && z_1(i)<0
                                      display case2
                                      
             
                                     k=i;
                                     m=k+1;
                                     
                                     
                                     for m=k+1:N
                                     t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     t_22=t_2_2;
                                     t_2_2=0;
                                     t_2(i)=t_22
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x_R(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                   
                                  elseif  z_1(i)*z_1(i+1)==0 && z_1(i)==0
                                      display case3
                                      k=i;
                                     m=k+1;
                                     
                                     
                                     for m=k+1:N
                                     t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     t_22=t_2_2;
                                     t_2_2=0;
                                     t_2(i)=t_22
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x_R(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                   
                                  elseif z_1(i)*z_1(i+1)==0 && z_1(i+1)==0
                                     display case4
                                     k=i-1;
                                     m=k+1;
                                     
                                     
                                         
                                     for m=k+1:N
                                     t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     
                                     
                                     t_21=t_2_1;
                                     t_2_1=0;
                                     t_2(i)=t_21
                                     Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                     T_x_R(:,:,i) = -(Q_x(:,:,i)*Q)/I
     
                                     end
     
                                     
         % t-2 calutaion ends                            
     
    else
        
        % t-1 and t-2 calutaion starts
        for v= 1:length(x_R)
           
         if x_R(:,v,i)>=0
         display case6  
         k=i-1;
         m=k+1;
         t_1(:,v,i) = (b/2)*n(i)*((x_R(:,v, i )-z_g).^2-(z(m))^2)
                                        
                                     
                                     
                                         
                                        for m=k+1:N
                                        t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                        end
                                     
                                     
                                        t_21=t_2_1
                                        t_2_1=0;
                                        t_2(i)=t_21
                                        Q_x(:,v,i) = t_1(:,v,i) + t_2(i); 
    
                                        T_x_R(:,v,i) = -(Q_x(:,v,i)*Q)/I
         
         else
         display case7
         k=i;
         m=k+1;
         t_1(:,v,i) = (b/2)*n(i)*((x_R(:,v, i )-z_g).^2-(z(m))^2)
                                        k=i;
                                        m=k+1;
                                     
                                     
                                        for m=k+1:N
                                        t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                        end
                                        t_22=t_2_2
                                        t_2_2=0;
                                        t_2(i)=t_22
                                        Q_x(:,v,i) = t_1(:,v,i) + t_2(i); 
    
                                        T_x_R(:,v,i) = -(Q_x(:,v,i)*Q)/I
       
         end
         
        v=v+1; 
        
        end
        
        
       % t-1and t-2  calutaion ends
    end  
     
                                     
                                   
end



M_R=max(T_x_R(:,:,:))

M_L=max(T_x_L(:,:,:))


% PLOTTING %   

% figure (1)
% plot (  T_x(:,:,1), x_R(:,:, 1) ) ;
 
% figure (1)
% plot (  T_x(:,:,1), x_R(:,:, 1) , T_x(:,:,2), x_R(:,:, 2) );
 
%  
 figure (1)
 plot (  T_x_R(:,:,1), x_R(:,:, 1) , T_x_R(:,:,2), x_R(:,:, 2) ,T_x_R(:,:,3) , x_R(:,:, 3),T_x_L(:,:,1), x_L(:,:, 1) , T_x_L(:,:,2), x_L(:,:, 2) ,T_x_L(:,:,3) , x_L(:,:, 3) );
  

% figure (2)
% plot (  T_x(:,:,1), x_R(:,:, 1) , T_x(:,:,2), x_R(:,:, 2) ,T_x(:,:,3) ,x_R(:,:, 3) , T_x(:,:,4) ,x_R(:,:, 4) ) ;

% figure (2)
% plot (  T_x(:,:,1), x_R(:,:, 1) , T_x(:,:,2), x_R(:,:, 2) ,T_x(:,:,3) ,x_R(:,:, 3) , T_x(:,:,4) ,x_R(:,:, 4) ) ;
 
% figure (2)
% plot (  T_x(:,:,1), x_R(:,:, 1) , T_x(:,:,2), x_R(:,:, 2) ,T_x(:,:,3) ,x_R(:,:, 3) , T_x(:,:,4) ,x_R(:,:, 4)  ,T_x(:,:,5) , x_R(:,:, 5) ) ;

% figure (1)
% plot (  T_x(:,:,1), x_R(:,:, 1) , T_x(:,:,2), x_R(:,:, 2) ,T_x(:,:,3) ,x_R(:,:, 3) , T_x(:,:,4) ,x_R(:,:, 4)  ,T_x(:,:,5) , x_R(:,:, 5) ,T_x(:,:,6) , x_R(:,:, 6) ) ;



T=[T_x_R(:,:,1) T_x_R(:,:,1) T_x_R(:,:,1)]

save('file.dat','T')













