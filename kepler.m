function [EH1,EH2,nu2]=kepler(t1,nu1,t2,a,e)

%Constants
mue=398600;
tolerance=1e-10;
nit=100;

%%%Elliptic Orbit%%%
if (e)<(1),
    %Inital Eccentric Anomaly
    E1=2*atan(tan(nu1/2)*((1-e)/(1+e))^0.5);
    if (E1)<(0)
      E1=2*pi+E1;
    end
    EH1=E1;
        
    %Final Eccentric Anomaly
    %N.I. Inital Guess
    E2_new=E1+(mue/a^3)^0.5*(t2-t1);
        for i=1:nit;
           E2_old=E2_new;
           fn=E2_old-E1-e*(sin(E2_old)-sin(E1))-(mue/a^3)^0.5*(t2-t1);
           fnd=1-e*cos(E2_old);
           E2_new=E2_old-fn/fnd;
           if (abs(E2_new-E2_old))<(tolerance)
                break
           end
        end
    E2=E2_new;
    no=round(E2/(2*pi)-0.5);
    E2=E2-no*2*pi;
    EH2=E2;

    %True Anomaly
    nu2=2*atan(tan(E2/2)*((1+e)/(1-e))^0.5);
    if (nu2)<(0);
        nu2=2*pi+nu2;
    end
end

%%%Hyperbolic Orbit%%%
if (e)>1,
    %Inital Hyperboliv Anomaly
    H1=2*atanh(tan(nu1/2)*((e-1)/(e+1))^0.5);
    EH1=H1;
        
    %Final Eccentric Anomaly
    %N.I. Inital Guess
    H2_new=H1;
        for i=1:nit;
           H2_old=H2_new;
           fn=-(H2_old-H1)+e*(sinh(H2_old)-sinh(H1))-(mue/a^3)^0.5*(t2-t1);
           fnd=1-e*cosh(H2_old);
           H2_new=H2_old-fn/fnd;
           if (abs(H2_new-H2_old))<(tolerance)
                break
           end
        end
    H2=H2_new;
    EH2=H2;

    %Final True Anomaly
    nu2=2*atan(tanh(H2/2)*((e+1)/(e-1))^0.5);
    if (nu2)<(0);
        nu2=2*pi+nu2;
    end
end
