function y = g(phi_ij,ae_zi,Rf,Rc)          

%{
Function:
    Determine whether the current cutter element is in cutting

Input:
    phi_ij: Radial angle
    ae_zi:  Cutting width
    Rf:     Cutting radius of the previous cutter element
    Rc:     Cutting radius of the current cutter element

Output:
    y = 1: cutting; y = 0: not cutting
%}

%% 
global up_or_down fe
if up_or_down == 1 % Up milling
    if Rf < Rc
        phi_st = -acos(Rf/Rc);
        phi_ex = pi/2 - asin( (Rf-ae_zi) / Rc );
    else
        phi_st = acos( (Rc^2 + fe^2 - Rf^2) / (2*Rc*fe) ) - pi/2;
        phi_ex = pi/2 - asin( (Rf-ae_zi) / Rc );
    end
    
elseif up_or_down == -1 % Down milling
    if Rf < Rc
        phi_st = pi/2 + asin( (Rf-ae_zi) / Rc );
        phi_ex = pi + acos(Rf/Rc);
    else
        phi_st = pi/2 + asin( (Rf-ae_zi) / Rc );
        phi_ex = 3*pi/2 - acos( (Rc^2 + fe^2 - Rf^2) / (2*Rc*fe) );
    end
end
phi_ij = mod(phi_ij+pi/2 , 2*pi) - pi/2;  % Adjust the radial angle to [-pi/2, 3pi/2)
if(phi_ij>=phi_st && phi_ij<=phi_ex)
    y=1;
else
    y=0;
end

end