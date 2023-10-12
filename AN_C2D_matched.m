function [D_z_result]=AN_C2D_matched(num,den,h,omega)

s_poles = roots(den); %Calculate the poles in the s-domain
s_zeroes = roots(num);%Calculate the zeroes in the s-domain

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

D_s_evaluate = RR_evaluate(Ds,1i*omega);  %Evaluate s = i*omega_bar for "prewarping"
D_z_evaluate = RR_evaluate(D_z, exp(1i*omega*h)); %Evaluate z = exp(i*omega_bar*h)

gain = D_s_evaluate/D_z_evaluate; %Scale D(z) to equal D(s)
D_z_result = gain*D_z;
end


