function P = P_func(x,y,x_v,y_v,rho_inf,p_inf,V_inf,phi)
    q_inf = 0.5*rho_inf*V_inf^2;
    Cp = 1 - (phi./V_inf).^2;
    
    P = p_inf + Cp.*q_inf;
end