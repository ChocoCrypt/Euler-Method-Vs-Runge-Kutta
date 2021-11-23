%------------------------------metodo de euler-------------------------------
%-----------------------------setup----------------------------------------
clear all
%defino el modelo SIR de infectados/recuperados
delta_t=1;
s=1 - 1e-6;
i=1e-6;
r=1-(s+i);
beta = 1.4247;
gamma = 0.14286;
total_pasos=60;
%inicio esto como una lista vacia
t=[0];
%iteracion metodo de euler
%resuelvo la ecuacion diferencial de 0  60, aunque en caso de que
for n=1:total_pasos
    t(n+1)=t(n)+delta_t;
    s(n+1)=s(n)-beta*s(n)*i(n)*delta_t;
    i(n+1)=i(n)+(beta*s(n)*i(n)-gamma*i(n))*delta_t;
    r(n+1)=r(n)+gamma*i(n)*delta_t;
end
%------------------grafico los valores-----------------------%

figure
plot(0:total_pasos,s)
hold on
plot(0:total_pasos,i)
hold on
plot(0:total_pasos,r)
legend("s" , "i" , "r")
%%
%-----------------------------setup----------------------------------------
hold on
clear all
%defino el modelo SIR de infectados/recuperados
delta_t=0.4;
s=1 - 1e-6;
i=1e-6;
r=1-(s+i);
beta = 1.4247;
gamma = 0.14286;
total_pasos=90;
%inicio esto como una lista vacia
t=[0];
%iteracion metodo de euler
%resuelvo la ecuacion diferencial de 0  90, aunque en caso de que
for n=1:total_pasos
    k_s1 = -beta*s(n)*i(n);
    k_i1 = (beta*s(n) * (i(n))) - (gamma*i(n));
    k_r1 = gamma*i(n);
    k_s2 = -beta* (s(n)+ delta_t*0.5 * k_s1) * (i(n) + k_i1*0.5 * delta_t);
    k_i2 = beta*(s(n) + 0.5 * delta_t * k_s1) *(i(n) + 0.5 * delta_t * k_i1) - gamma*(i(n) + 0.5*k_i1);
    k_r2 = gamma * (i(n) + 0.5 *delta_t * k_i1 );
    k_s3 = -beta*(s(n) + 0.5 * delta_t * k_s2 )*(i(n) + 0.5 * delta_t * k_i2);
    k_i3 = beta * (s(n) + 0.5 * delta_t * k_s2 )*(i(n)+0.5  * delta_t * k_i2) - gamma*(i(n)+0.5 * delta_t * k_i2);
    k_r3 = gamma * (i(n)+0.5 * delta_t * k_i2) ;
    k_s4 = -beta * (s(n) + delta_t * k_s3) * (i(n) + delta_t *k_i3); %BRE
    k_i4 = beta*(s(n) + delta_t * k_s3)*(i(n) + delta_t * k_i3) - gamma * (i(n) + delta_t * k_i3);
    k_r4 = gamma * (i(n) + delta_t * k_i3);
    s(n+1) = s(n) + delta_t * (k_s1 +2*k_s2 + 2*k_s3 + k_s4)/6;
    i(n+1) = i(n) + delta_t * (k_i1 +2*k_i2 + 2*k_i3 + k_i4)/6;
    r(n+1) = r(n) + delta_t * (k_r1 +2*k_r2 + 2*k_r3 + k_r4)/6;
end
figure
plot(0:total_pasos,s)
hold on
plot(0:total_pasos,i)
hold on
plot(0:total_pasos,r)
legend("s" , "i" , "r")