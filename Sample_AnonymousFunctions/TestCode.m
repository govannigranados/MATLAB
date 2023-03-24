%This Script will run tests on this set exercise set.
%===================================================
%Testing with scripts will verify the base code works before you make it
%into a reusable function.  You should be able to copy paste most of the
%work from the test script to the final product.

%============================
%First we will set some grids
%============================
Grid1 = 0:.1:1;
Grid2 = 0:.1:2;

%Run the tests
%Here we call our Sample grid functions
%----------------------------------------------------------------------
%Using @(x)FunctionName(const,x) as input forces MATLAB to evaluate the function
%object always using const as the first input and allowing x to varying.
%----------------------------------------------------------------------
%The function Sample_AnonymousFunctionGrid will now calculate the value
%on the given grid changing the arguement i over the loop.
for i=1:3
    Poly1(i,:)=Sample_AnonymousFunctionGrid(@(x)Input_FirstDegreePoly(i,x),Grid1);
    Poly2(i,:)=Sample_AnonymousFunctionGrid(@(x)Input_FirstDegreePoly(i,x),Grid2);
end
%=============
%Plot Results
%==============
%We wil use both [0,1] and [0,2] and verify by eye.
%[0,1]
figure()
hold on
for i =1:3
    plot(Grid1, Poly1(i,:),"LineWidth",2.0,"MarkerSize",20)
end
legend("x","2x","3x")
hold off
title("f(x)=cx")
%-------------------------
%[0,2]
figure()
hold on
for i =1:3
    plot(Grid2, Poly2(i,:),"LineWidth",2.0,"MarkerSize",20)
end
legend("x","2x","3x")
hold off
title("f(x)=cx")