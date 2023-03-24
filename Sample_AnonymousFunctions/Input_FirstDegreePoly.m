function Result = Input_FirstDegreePoly(Const,X)
%This is intented to be used as an anonymous function for:
%"Sample_AnonymousFunctionGrid"
%It should be used as @(x)Input_FirstDegreePoly(-,x)
%You can hardcode '-' when you call the anonymous function
%===========================================================
Result = Const*X;
end

