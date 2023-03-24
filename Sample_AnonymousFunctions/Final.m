function Final(F,Grid,Title)
%This Turns the script into a callable function.
%===================================================
%F is the function you are evaluating, call as 
% 1: @(i,x)function(i,x) if written yourself
% 2: @(i,x)handle(args) if matlab recognizes it (sine, exp etc)
%--------------------------------------------------------------------------
%Grid is a set of points you can evaluate over.  [0,1] etc in discrete form
%--------------------------------------------------------------------------
%Title is what you want your plot title to be, a string "f(x)=foo".
%===================================================
%Here we call our Sample grid functions
%Try: Final(@(i,x)Input_FirstDegreePoly(i,x), 0:.1:1,"f(x)=ix")
%Try: Final(@(i,x)sin(i*x), 0:.1:1, "f(x)=sin(ix)")
%Try: Final(@(i,x)exp(i*x), 0:.1:1,"f(x)=exp(ix)")
%We will add some other functions to this folder and try some more
%----------------------------------------------------------------------
for i=1:3
    %This trick will force the F function to do the following
    %Force i to be fixed
    % allow x to be evaluated as a variable
    Poly(i,:)=Sample_AnonymousFunctionGrid(@(x)F(i,x),Grid);
end
%Plot Results
%==============
%[0,1]
figure()
hold on
for i =1:3
    plot(Grid, Poly(i,:),"LineWidth",2.0,"MarkerSize",20)
end
legend("x","2x","3x")
hold off
title(Title)


