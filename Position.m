function[l3,L3,o3,o4,s] = Position(l1,l2,L,l4,l5,o2)

    X = [500;20/180*pi];
    tol = 1e-9;
    change = [inf;inf];
    
    l3 = X(1);
    o3  = X(2);
    L3 = L-l3;
    iter = 0;
    
    max_iter = 100;
    while(abs(change(1))>tol && abs(change(2))>tol)
        l3 = X(1);
        o3 = X(2);
        J = [-cos(o3),l3*sin(o3);-sin(o3),-l3*cos(o3)];
        F = [l1+l2*cos(o2)-l3*cos(o3);l2*sin(o2)-l3*sin(o3)];
        if(abs(det(J))<1e-6 )
            fprintf("Approaching singular jacobian")
            break;
        end
        change = -1*pinv(J)*F;
        X = X +change;
        iter = iter + 1;
        if(iter>=max_iter)
            fprintf("overflow")
            break;
        end
        
    end
    l3 = X(1);
    o3  = X(2);
    L3 = L-l3;

    Y = [70;60/180*pi];
    tol = 1e-9;
    change = [inf;inf];
    s = Y(1);
    o4 = Y(2);
    iter = 0;
    max_iter = 100;
    while(abs(change(1))>tol && abs(change(2))>tol)
        s = Y(1);
        o4 = Y(2);
        J = [0,l4*sin(o4);1,-l4*cos(o4)];
        F = [l2*cos(o2) + L3*cos(o3)-l4*cos(o4)-l5;l2*sin(o2) + L3*sin(o3)-l4*sin(o4)+s];
        if(abs(det(J))<1e-6 )
            fprintf("Approaching singular jacobian")
            break;
        end
        change = -1*pinv(J)*F;
        Y = Y +change;
        iter = iter + 1;
        if(iter>=max_iter)
            fprintf("overflow")
            break;
        end
    
    end

    s = Y(1);
    o4 = Y(2);

end

