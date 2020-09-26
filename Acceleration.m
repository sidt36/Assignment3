function[s_2,alpha_4] = Acceleration(l1,l2,L,l4,l5,o2,w2)
[l3,L3,o3,o4,~] = Position(l1,l2,L,l4,l5,o2);
[l3_1,L3_1,w4,~,w3] = Velocity(l1,l2,L,l4,l5,o2,w2);
J = [-cos(o3),l3*sin(o3);-sin(o3),-l3*cos(o3)];


C = w2*w2*[l2*cos(o2);l2*sin(o2)] + 2*l3_1*w3*[-sin(o3);cos(o3)] + w3*w3*[-l3*cos(o3);-l3*sin(o3)];
X2 = pinv(J)*C;

alpha_3 = X2(2);
l3_2 = X2(1);
L3_2 = -l3_2;


C1 = w2^2 * l2 * [cos(o2);sin(o2)] + L3_1*w3*[2*sin(o3);-2*cos(o3)] + L3_2*[-cos(o3);-sin(o3)] + alpha_3*L3*[sin(o3);-cos(o3)]+ L3*(w3)^2*[cos(o3);sin(o3)] + l4*(w4)^2*[-cos(o4);-sin(o4)];
J = [0,l4*sin(o4);1,-l4*cos(o4)];
Y2 =  pinv(J) * C1;




s_2 = Y2(1);
alpha_4 = Y2(2);
end