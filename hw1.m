%% Problem 1
disp('https://github.com/Arjun-Naageshwaran')   %Profile with all my repos and pinned projects
disp('https://github.com/Arjun-Naageshwaran/UCSD-MAE-144.git')   %Branch of UCSD-Classes repo which contains my MAE 144 assignments

%% Problem 2a
b = RR_poly([-2 2 -5 5], 1);   %Numerator of the plant
a = RR_poly([-1 1 -3 3 -6 6], 1);   %Denominator of the plant
f_1 = RR_poly([-1 -1 -3 -3 -6 -6], 1);   %Target poles

[x_1,y_1] = RR_Diophantine(a,b,f_1)   %Return smallest order integer solution
test_1 = trim(a*x_1+b*y_1)
residual_1 = norm(f_1-test_1)

%% Problem 2b
b = RR_poly([-2 2 -5 5], 1);   %Numerator of the plant
a = RR_poly([-1 1 -3 3 -6 6], 1);   %Denominator of the plant
f_2 = RR_poly([-1 -1 -3 -3 -6 -6 -20 -20 -20 -20 -20 -20], 1);   %Target poles

[x_2,y_2] = RR_Diophantine(a,b,f_2)
test_2 = trim(a*x_2+b*y_2)
residual_2 = norm(f_2-test_2)
%The controller in 2a is improper. This is because y_1 has a degree of 5
%which exceeds that of x_1, 3. I changed the target function by adding k=6
%additional poles which creates a proper controller. Now y_2 with a degree
%of 5 does not exceed x_2 which has a degree of 6.

%% Problem 3
%AN_C2D_matched
s_poles = roots(den); %Calculate the poles in the s-domain
s_zeroes = roots(num);%Calculate the zeroes in the s-domain
%omega_bar = 0   %Arbitrary value that can be changed between bounds 0 < omega_bar < pi/h

z_poles = exp(s_poles*h); %Use equation z=e^(s*h) to map D(s) poles to D(z)
z_zeroes = exp(s_zeroes*h); %Use equation z=e^(s*h) to map D(s) zeroes to D(z)
l = length(z_poles)-length(z_zeroes);

%Add infinite zeroes by substracting #poles - #zeroes and then mapping them
%into a matrix
for i = 1:l 
    z_zeroes(end+1,1) = (-1);
end

l = length(z_poles)-length(z_zeroes);

%Condition of strict causality
if l > 0
    if isinf(z_poles) 
        z_zeroes(end) = inf;  %Map 1 zero to z = infity in D(z)
    end
%Condition of semi-causality
elseif l == 0 
    z_zeroes(end) = []; %Take away one zero to make it strictly causal
else
end

s_numerator = poly(s_zeroes); 
s_denominator = poly(s_poles);
Ds = RR_tf(s_numerator,s_denominator);

z_numerator = poly(z_zeroes);
z_denominator = poly(z_poles);
D_z = RR_tf(z_numerator,z_denominator);

D_s_evaluate = RR_evaluate(Ds,1i*omega_bar);  %Evaluate s = i*omega_bar for "prewarping"
D_z_evaluate = RR_evaluate(D_z, exp(1i*omega_bar*h)); %Evaluate z = exp(i*omega_bar*h)

gain = D_s_evaluate/D_z_evaluate; %Scale D(z) to equal D(s)
D_z_result = gain*D_z;
