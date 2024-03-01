function [m,c,k] = modal(f,damping,u)

%{
Function:
    Calculate modal parameters

Input:
    f:          Damped modal frequency
    damping:    Damping ratio
    u:          Imaginary value of the residue

Output:
    m:	Modal mass
    c: 	Modal damping
    k:	Modal stiffness
%}

%%
u = abs(u);
wd = f*2*pi;                    % Damped modal angular frequency
wn = wd./sqrt(1-damping.^2);    % Undamped modal angular frequency
m = 1./(2*u.*wn);     
k = wn.^2.*m;                
c = 2*damping.*sqrt(m.*k); 

end