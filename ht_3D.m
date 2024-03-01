function ht = ht_3D(z_Num)

%{
Function:
    Calculate h(t)

Input:
    z_Num: Number of tool discretizations in the axial direction

Output:
    ht: Matrix h(t)
%}

%% Calculate hxx, hxy, hyx, hyy
global N Kt Kr Ka mm dz DD r beta ae lambda rou
ht = zeros(3,3,mm+1);
[a11,a12,a13,a21,a22,a23,a31,a32,a33] = deal(zeros(mm+1,1));
for i = 1: mm+1             % Traverse the spindle rotation period
    dtr = 2*pi/mm;          % Tool rotation angle
    for j=1:N               % Traverse the tool teeth
        for iz = 1:z_Num    % Traverse the axial height of the tool
            zi = dz*iz;
            phi_ij = i*dtr+(j-1)*2*pi/N-2*zi*tan(beta)/DD;      % Radial angle of the cutter element
            
            if zi < r
                k = acos(1-zi/r);                               % Axial contact angle of the cutter element
                R_zi = DD/2-r+sqrt(r^2-(r-zi)^2);               % Cutting radius of the cutter element
            else
                k = pi/2;
                R_zi = DD/2;
            end
           %% 
            ae_zi = ae; % Cutting width of the cutter element
            
            phi_zi_j_runout = lambda + 2*zi*tan(beta)/DD - (j-1)*2*pi/N;    % Angle between the current cutter element and the runout direction
            Rc = sqrt( R_zi^2 + rou^2 - 2*R_zi*rou*cos(phi_zi_j_runout) );  % Cutting radius of the current cutter element
            phi_zi_jf_runout = lambda + 2*zi*tan(beta)/DD - (j)*2*pi/N;     % Angle between the previous cutter element and the runout direction
            Rf= sqrt( R_zi^2 + rou^2 - 2*R_zi*rou*cos(phi_zi_jf_runout) );  % Cutting radius of the previous cutter element
            
            g_ij = g(phi_ij,ae_zi,Rf,Rc); % Cutting judgment function of the cutter element
            sp = sin(phi_ij); cp = cos(phi_ij); sk = sin(k); ck = cos(k);
            a11(i) = a11(i) + g_ij*(- Kt*sk*sp*cp   - Kr*sk^2*sp^2      - Ka*sk*ck*sp^2)    *dz/sk;
            a12(i) = a12(i) + g_ij*(- Kt*sk*cp^2    - Kr*sk^2*sp*cp     - Ka*sk*ck*sp*cp)   *dz/sk;
            a13(i) = a13(i) + g_ij*(+ Kt*ck*cp      + Kr*sk*ck*sp       + Ka*ck^2*sp)       *dz/sk;
            
            a21(i) = a21(i) + g_ij*(+ Kt*sk*sp^2    - Kr*sk^2*sp*cp     - Ka*sk*ck*sp*cp)   *dz/sk;
            a22(i) = a22(i) + g_ij*(+ Kt*sk*sp*cp   - Kr*sk^2*cp^2      - Ka*sk*ck*cp^2)    *dz/sk;
            a23(i) = a23(i) + g_ij*(- Kt*ck*sp      + Kr*sk*ck*sp       + Ka*ck^2*cp)       *dz/sk;
            
            a31(i) = a31(i) + g_ij*(                + Kr*sk*ck*sp       - Ka*sk^2*sp)       *dz/sk;
            a32(i) = a32(i) + g_ij*(                + Kr*sk*ck*cp       - Ka*sk^2*cp)       *dz/sk;
            a33(i) = a33(i) + g_ij*(                - Kr*ck^2           + Ka*sk*ck)         *dz/sk;
        end
    end
    ht(:,:,i) = [a11(i),    a12(i),     a13(i)
                 a21(i),    a22(i),     a23(i)
                 a31(i),    a32(i),     a33(i)];

end
end