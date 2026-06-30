%% Trabajo I - Econometría II
% Estudiantes: Fernanda Anguita, Santiago García, Sebastian Mejías, Sergio
% Salazar y Leonardo Siles

clc, clear;
cd('/Users/sergiosalazar/Documents/Magister/Segundo Semestre/Econometría II/Tareas/hmw1_24');

%% Cargar la data
data=readmatrix('hmw1_24.xlsx');

y1=data(:,1);
y2=data(:,2);
y3=data(:,3);

%% Determinación de (p,q) por Hannan-Rissanen
L=25; %Maximo valor al que se recorren las matrices de rezagos
g=L; %Valor arbitrario del orden del proceso auxiliar

%Creando los rezagos de cada serie (L rezagos)
rez1=zeros(size(y1,1),L);

for j=1:L
    for i=1:(size(y1,1)-j)
        rez1(i+j,j)=y1(i,:);
    end
end

% Paso (1): Estimación auxiliar de un AR(j) 
autoe1=zeros(length(y1),L);

for j=1:g %OLS con g suficientemente grande
    [~,emat1]=OLS(y1,rez1(:,1:j)); %Cuidar que g=L por lo menos
    e1(:,j)=emat1; %Almacena los residuos de cada g regresión
end

for i=1:L %Creamos matrices de rezagos del residuo
      for j=1:(length(y1)-i)
        autoe1(i+j,i)=e1(j,g);
    end
end

%Paso (2): Se guardan los residuos del AR(j) para testearlos H0:Ruido
%blanco
%Se construyen autocorrelaciones con un valor s arbitrario que es el valor
%del rezago con el que se correlaciona el residuo

T=length(y1); %Tamaño de muestra
alpha=0.05; %Valor de significancia
s=L;%Máximo hasta el valor de g

for i=0:L-1
    den1(i+1)=sum(emat1.*autoe1(:,i+1)); %Denominador
    num1=sum(emat1.^2);
   rs1(i+1)=den1(i+1)/num1; %Autocorrelación serial de residuos. Un valor por cada i
    Q1(i+1)=(T*(T+2)*(sum(rs1(:,i+1:end).^2))/(T-i+1)); %Test de Lljung-Box
    vc1(i+1)=chi2inv(1-alpha,i+1); 
end

    if Q1(g)>vc1(g)
       disp('No se acepta H0, entonces los residuos no son ruido blanco');
    else 
       disp('Se acepta H0, entonces los residuos son ruido blanco');
    end
 % El mensaje que sale toma en cuenta el valor más alto 
%pval=1-chi2cdf(Q1,s);

%Paso (3): Se guardan los valores del rezago resultante como white noise
ma1=[emat1 autoe1]; %Admite ceros porque probamos modelos puros
ar1=[rez1]; %Admite ceros porque probamos modelos puros
vec1=[ones(length(y1),1) ar1 ma1]; %vec1=[1 0 ar(L) 0 ma(L) ]
maxlag=L; %El número de esto depende el núm de modelos a testear 

for p=0:maxlag
    for q=0:maxlag
         if p == 0 && q == 0
            X1 = vec1(:, 1); % Solo el término constante
         elseif q == 0
            X1 = vec1(:, [1, 2:p+1]); % Solo AR terms
         elseif p==0
             X1=vec1(:,[1, maxlag+3:maxlag+2+q]); %Solo MA terms
        else
            X1 = vec1(:, [1, 2:p+1, maxlag+3:maxlag+2+q]); % AR and MA terms
         end
        [beta,res1,ss,kk]=OLS(y1(maxlag+1:end,:),X1(maxlag+1:end,:));
        phi1{p+1,q+1}=beta;
        emodel1{p+1,q+1}=res1;
        s2_1{p+1,q+1}=ss;
        K1{p+1,q+1}=kk; 
    end
end 
for p=1:maxlag+1
    for q=1:maxlag+1
         bic1(p,q)=log(s2_1{p,q})+(K1{p,q}/(T-maxlag))*log(T-maxlag);
        hqc1(p,q)=log(s2_1{p,q})+2*(K1{p,q}/(T-maxlag))*log(log(T-maxlag));
    end
end

[min_bic1, idx_bic1] = min(bic1(:));
[row_bic1, col_bic1] = ind2sub(size(bic1), idx_bic1);
[min_hqc1, idx_hqc1] = min(hqc1(:));
[row_hqc1, col_hqc1] = ind2sub(size(hqc1), idx_hqc1);
fprintf('Modelo 1 con BIC mínimo: p = %d, q = %d, BIC = %.4f\n', row_bic1, col_bic1, min_bic1);
fprintf('Modelo 1 con HQC mínimo: p = %d, q = %d, HQC = %.4f\n', row_hqc1, col_hqc1, min_hqc1);
if row_bic1 == row_hqc1 && col_bic1 == col_hqc1
    fprintf('El modelo  1 con el mejor valor tanto en BIC como en HQC es: p = %d, q = %d\n', row_bic1, col_bic1);
else
    fprintf('El mejor modelo 1 según BIC es: p = %d, q = %d\n', row_bic1, col_bic1);
    fprintf('El mejor modelo 1 según HQC es: p = %d, q = %d\n', row_hqc1, col_hqc1);
end %(LIMPIO)

%% Aplicamos el mismo algoritmo para la serie 2
clc;
L=25; %Maximo valor al que se recorren las matrices de rezagos
g=L; %Valor arbitrario del orden del proceso auxiliar

%Creando los rezagos de cada serie (L rezagos)
rez2=zeros(size(y2,1),L);

for j=1:L
    for i=1:(size(y2,1)-j)
        rez2(i+j,j)=y2(i,:);
    end
end

% Paso (1): Estimación auxiliar de un AR(j) 
autoe2=zeros(length(y2),L);

for j=1:g %OLS con g suficientemente grande
    [~,emat2]=OLS(y2,rez2(:,1:j)); %Cuidar que g=L por lo menos
    e2(:,j)=emat2; %Almacena los residuos de cada g regresión
end

for i=1:L %Creamos matrices de rezagos del residuo
      for j=1:(length(y2)-i)
        autoe2(i+j,i)=e2(j,g);
    end
end

%Paso (2): Se guardan los residuos del AR(j) para testearlos H0:Ruido
%blanco
%Se construyen autocorrelaciones con un valor s arbitrario que es el valor
%del rezago con el que se correlaciona el residuo

T=length(y2); %Tamaño de muestra
alpha=0.05; %Valor de significancia
s=L;%Máximo hasta el valor de g

for i=0:L-1
    den2(i+1)=sum(emat2.*autoe2(:,i+1)); %Denominador
    num2=sum(emat2.^2);
   rs2(i+1)=den2(i+1)/num2; %Autocorrelación serial de residuos. Un valor por cada i
    Q2(i+1)=(T*(T+2)*(sum(rs2(:,i+1:end).^2))/(T-i+1)); %Test de Lljung-Box
    vc2(i+1)=chi2inv(1-alpha,i+1); 
end

    if Q2(g)>vc2(g)
       disp('No se acepta H0, entonces los residuos no son ruido blanco');
    else 
       disp('Se acepta H0, entonces los residuos son ruido blanco');
    end
 % El mensaje que sale toma en cuenta el valor más alto 
%pval=1-chi2cdf(Q1,s);

%Paso (3): Se guardan los valores del rezago resultante como white noise
ma2=[emat2 autoe2]; %Admite ceros porque probamos modelos puros
ar2=[rez2]; %Admite ceros porque probamos modelos puros
vec2=[ones(length(y2),1) ar2 ma2]; %vec1=[1 0 ar(L) 0 ma(L) ]
maxlag=L; %El número de esto depende el núm de modelos a testear 

for p=0:maxlag
    for q=0:maxlag
         if p == 0 && q == 0
            X2 = vec2(:, 1); % Solo el término constante
         elseif q == 0
            X2 = vec2(:, [1, 2:p+1]); % Solo AR terms
         elseif p==0
             X2=vec2(:,[1, maxlag+3:maxlag+2+q]); %Solo MA terms
        else
            X2 = vec2(:, [1, 2:p+1, maxlag+3:maxlag+2+q]); % AR and MA terms
         end
        [beta2,res2,ss2,kk2]=OLS(y2(maxlag+1:end,:),X2(maxlag+1:end,:));
        phi2{p+1,q+1}=beta2;
        emodel2{p+1,q+1}=res2;
        s2_2{p+1,q+1}=ss2;
        K2{p+1,q+1}=kk2; 
    end
end 
for p=1:maxlag+1
    for q=1:maxlag+1
         bic2(p,q)=log(s2_2{p,q})+(K2{p,q}/(T-maxlag))*log(T-maxlag);
        hqc2(p,q)=log(s2_2{p,q})+2*(K2{p,q}/(T-maxlag))*log(log(T-maxlag));
    end
end

[min_bic2, idx_bic2] = min(bic2(:));
[row_bic2, col_bic2] = ind2sub(size(bic2), idx_bic2);
[min_hqc2, idx_hqc2] = min(hqc2(:));
[row_hqc2, col_hqc2] = ind2sub(size(hqc2), idx_hqc2);
fprintf('Modelo 1 con BIC mínimo: p = %d, q = %d, BIC = %.4f\n', row_bic2, col_bic2, min_bic2);
fprintf('Modelo 1 con HQC mínimo: p = %d, q = %d, HQC = %.4f\n', row_hqc2, col_hqc2, min_hqc2);
if row_bic2 == row_hqc2 && col_bic2 == col_hqc2
    fprintf('El modelo  1 con el mejor valor tanto en BIC como en HQC es: p = %d, q = %d\n', row_bic2, col_bic2);
else
    fprintf('El mejor modelo 1 según BIC es: p = %d, q = %d\n', row_bic2, col_bic2);
    fprintf('El mejor modelo 1 según HQC es: p = %d, q = %d\n', row_hqc2, col_hqc2);
end %(LIMPIO)

%% Aplicamos el mismo algoritmo para la serie 3
clc;
L=25; %Maximo valor al que se recorren las matrices de rezagos
g=L; %Valor arbitrario del orden del proceso auxiliar

%Creando los rezagos de cada serie (L rezagos)
rez3=zeros(size(y3,1),L);

for j=1:L
    for i=1:(size(y3,1)-j)
        rez3(i+j,j)=y3(i,:);
    end
end

% Paso (1): Estimación auxiliar de un AR(j) 
autoe3=zeros(length(y3),L);

for j=1:g %OLS con g suficientemente grande
    [~,emat3]=OLS(y3,rez3(:,1:j)); %Cuidar que g=L por lo menos
    e3(:,j)=emat3; %Almacena los residuos de cada g regresión
end

for i=1:L %Creamos matrices de rezagos del residuo
      for j=1:(length(y3)-i)
        autoe3(i+j,i)=e3(j,g);
    end
end

%Paso (2): Se guardan los residuos del AR(j) para testearlos H0:Ruido
%blanco
%Se construyen autocorrelaciones con un valor s arbitrario que es el valor
%del rezago con el que se correlaciona el residuo

T=length(y3); %Tamaño de muestra
alpha=0.05; %Valor de significancia
s=L;%Máximo hasta el valor de g

for i=0:L-1
    den3(i+1)=sum(emat3.*autoe3(:,i+1)); %Denominador
    num3=sum(emat3.^2);
   rs3(i+1)=den3(i+1)/num3; %Autocorrelación serial de residuos. Un valor por cada i
    Q3(i+1)=(T*(T+2)*(sum(rs3(:,i+1:end).^2))/(T-i+1)); %Test de Lljung-Box
    vc3(i+1)=chi2inv(1-alpha,i+1); 
end

    if Q3(g)>vc3(g)
       disp('No se acepta H0, entonces los residuos no son ruido blanco');
    else 
       disp('Se acepta H0, entonces los residuos son ruido blanco');
    end
 % El mensaje que sale toma en cuenta el valor más alto 
%pval=1-chi2cdf(Q1,s);

%Paso (3): Se guardan los valores del rezago resultante como white noise
ma3=[emat3 autoe3]; %Admite ceros porque probamos modelos puros
ar3=[rez3]; %Admite ceros porque probamos modelos puros
vec3=[ones(length(y3),1) ar3 ma3]; %vec1=[1 0 ar(L) 0 ma(L) ]
maxlag=L; %El número de esto depende el núm de modelos a testear 

for p=0:maxlag
    for q=0:maxlag
         if p == 0 && q == 0
            X3 = vec3(:, 1); % Solo el término constante
         elseif q == 0
            X3 = vec3(:, [1, 2:p+1]); % Solo AR terms
         elseif p==0
             X3=vec3(:,[1, maxlag+3:maxlag+2+q]); %Solo MA terms
        else
            X3 = vec3(:, [1, 2:p+1, maxlag+3:maxlag+2+q]); % AR and MA terms
         end
        [beta3,res3,ss3,kk3]=OLS(y3(maxlag+1:end,:),X3(maxlag+1:end,:));
        phi3{p+1,q+1}=beta3;
        emodel3{p+1,q+1}=res3;
        s2_3{p+1,q+1}=ss3;
        K3{p+1,q+1}=kk3; 
    end
end 
for p=1:maxlag+1
    for q=1:maxlag+1
         bic3(p,q)=log(s2_3{p,q})+(K3{p,q}/(T-maxlag))*log(T-maxlag);
        hqc3(p,q)=log(s2_3{p,q})+2*(K3{p,q}/(T-maxlag))*log(log(T-maxlag));
    end
end

[min_bic3, idx_bic3] = min(bic3(:));
[row_bic3, col_bic3] = ind2sub(size(bic3), idx_bic3);
[min_hqc3, idx_hqc3] = min(hqc3(:));
[row_hqc3, col_hqc3] = ind2sub(size(hqc3), idx_hqc3);
fprintf('Modelo 1 con BIC mínimo: p = %d, q = %d, BIC = %.4f\n', row_bic3, col_bic3, min_bic3);
fprintf('Modelo 1 con HQC mínimo: p = %d, q = %d, HQC = %.4f\n', row_hqc3, col_hqc3, min_hqc3);
if row_bic3 == row_hqc3 && col_bic3 == col_hqc3
    fprintf('El modelo  1 con el mejor valor tanto en BIC como en HQC es: p = %d, q = %d\n', row_bic3, col_bic3);
else
    fprintf('El mejor modelo 1 según BIC es: p = %d, q = %d\n', row_bic3, col_bic3);
    fprintf('El mejor modelo 1 según HQC es: p = %d, q = %d\n', row_hqc3, col_hqc3);
end %(LIMPIO)