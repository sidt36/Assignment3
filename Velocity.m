function[l3_1,L3_1,w4,s_1,w3] = Velocity(l1,l2,L,l4,l5,o2,w2)

[l3,L3,o3,o4,s] = Position(l1,l2,L,l4,l5,o2);

J = [-cos(o3),l3*sin(o3);-sin(o3),-l3*cos(o3)];
w2 = -7*pi;
B = w2*[l2*sin(o2);-l2*cos(o2)];
X1 = pinv(J) * B;

w3 = X1(2);
l3_1 = X1(1);
L3_1 = -l3_1;


B1 = w3 *L3 *[sin(o3);-cos(o3)] + l2*w2*[sin(o2);-cos(o2)] + L3_1*[-cos(o3);-sin(o3)];
J = [0,l4*sin(o4);1,-l4*cos(o4)];
Y1 =  pinv(J) * B1;
Y1;
w4 = Y1(2);
s_1 = Y1(1);

end