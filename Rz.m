function R = Rz(theta)

%{
Function:
    Calculate the rotation matrix around Z axis

Input:
    theta: Rotation angle

Output:
    R: Rotation matrix
%}

%%
R = [cos(theta), -sin(theta), 0
     sin(theta), cos(theta) , 0
     0,          0          , 1];
end