%% Problem 1
disp('https://github.com/Arjun-Naageshwaran')   %Profile with all my repos and pinned projects
disp('https://github.com/Arjun-Naageshwaran/UCSD-Classes-Arjun.git')   %Branch of UCSD-Classes repo which contains my MAE 144 assignments

%% Problem 2a
a = RR_poly([-1 1 -3 3 -6 6], 1);
b = RR_poly([-2 2 -5 5], 1);
f_1 = RR_poly([-1 -1 -3 -3 -6 -6], 1);

[x_1,y_1] = RR_Diophantine(a,b,f_1)
test_1 = trim(a*x_1+b*y_1)
residual_1 = norm(f_1-test_1)

%% Problem 2b
a = RR_poly([-1 1 -3 3 -6 6], 1);
b = RR_poly([-2 2 -5 5], 1);
f_2 = RR_poly([-1 -1 -3 -3 -6 -6 -20 -20 -20 -20 -20 -20], 1);

[x_2,y_2] = RR_Diophantine(a,b,f_2)
test_2 = trim(a*x_2+b*y_2)
residual_2 = norm(f_2-test_2)
%The controller in 2a is improper. This is because y_1 has a degree of 5
%which exceeds that of x_1, 3. I changed the target function by adding k=6
%additional poles which creates a proper controller. Now y_2 with a degree
%of 5 does not exceed x_2 which has a degree of 6.