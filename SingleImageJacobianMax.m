function [NKMATJACOB, NFRAMEJACOB] = SingleImageJacobianMax(KMatrix ,NANGLE,... 
         NTRANSLATION, NCORRESPOND, NCONSENSUS)

    
%Counting the number of terms in NCONSENSUS
n = 0;
for j = 1 : length(NCONSENSUS)
    if NCONSENSUS(j) >= 1
       n = n + 1;
    end
end

%Generate a NKMATJACOB, a 2n x 5 matrix
%Generate a NFRAMEJACOB, a 2n x 6 matrix
NKMATJACOB = zeros(2*n, 5);
NFRAMEJACOB= zeros(2*n, 6);

perturbaion = 1.001;
for a = 1 : n
    %Extract corresponding point from NCORRESPOND
    b = NCONSENSUS(a);
    XY = NCORRESPOND(3:4, b);
    UV = NCORRESPOND(1:2, b);
    
    %Calculate perspectivies by finding rotaional matrix from Angle-axis
    ROTANG = norm(NANGLE);
    ROTAXIS = NANGLE/norm(NANGLE);
    R = RodriguesRotation(ROTAXIS, ROTANG);
    P = [R(1:3,1:2) NTRANSLATION]; %Perspectivities
    
    UVest = KMatrix * P * [XY;1];
    UVest = UVest(1:2)/UVest(3);
    
    k1 = KMatrix;
    k1(1,1) = KMatrix(1,1)*perturbaion;
    k2 = KMatrix;
    k2(1,2) = KMatrix(1,2)*perturbaion;
    k3 = KMatrix;
    k3(1,3) = KMatrix(1,3)*perturbaion;
    k4 = KMatrix; 
    k4(2,2) = KMatrix(2,2)*perturbaion;
    k5 = KMatrix;
    k5(2,3) = KMatrix(2,3)*perturbaion;
    
    fk1 = k1 * P *[XY; 1];
    fk2 = k2 * P *[XY; 1];
    fk3 = k3 * P *[XY; 1];
    fk4 = k4 * P *[XY; 1];
    fk5 = k5 * P *[XY; 1];
    
    dk1 = (fk1(1:2)/fk1(3)-UVest)/(KMatrix(1,1)*(perturbaion-1));
    dk2 = (fk2(1:2)/fk2(3)-UVest)/(KMatrix(1,2)*(perturbaion-1));
    dk3 = (fk3(1:2)/fk3(3)-UVest)/(KMatrix(1,3)*(perturbaion-1));
    dk4 = (fk4(1:2)/fk4(3)-UVest)/(KMatrix(2,2)*(perturbaion-1));
    dk5 = (fk5(1:2)/fk5(3)-UVest)/(KMatrix(2,3)*(perturbaion-1));
    
    NKMATJACOB(2*a-1:2*a,:) = [ dk1 dk2 dk3 dk4 dk5 ];
    
    for g = 1:3
        A = NANGLE;
        A(g) = NANGLE(g)*perturbaion;
        
        ROTANG = norm(A);
        ROTAXIS = A/norm(A);
        R = RodriguesRotation(ROTAXIS, ROTANG);
        P = [R(1:3,1:2) NTRANSLATION]; %Perspectivities
        
        fk = KMatrix * P* [XY; 1];
        
        dA = (fk(1:2)/fk(3) - UVest)/(NANGLE(g)*(perturbaion-1));
        
        NFRAMEJACOB( 2*a-1 : 2*a, g ) = dA;
    end
    
    for t = 1:3
        T = NTRANSLATION;
        T(t) = NTRANSLATION(t)*perturbaion;
        
        ROTANG = norm(NANGLE);
        ROTAXIS = NANGLE/norm(NANGLE);
        R = RodriguesRotation(ROTAXIS, ROTANG);
        P = [R(1:3,1:2) T]; %Perspectivities
        
        fk = KMatrix * P* [XY; 1];
        
        dT = (fk(1:2)/fk(3) - UVest)/ (NTRANSLATION(t)*(perturbaion-1));
        
        NFRAMEJACOB( 2*a-1 : 2*a, t+3 ) = dT;
    end
    
end
end
        
        
        
        
        
    
    
    
    
    
    
    
    
    