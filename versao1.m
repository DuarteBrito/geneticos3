% perfil inicial
perfilb=importdata('naca0009.txt', ' ');
x=perfilb(:,1);
y=perfilb(:,2);
extra = perfilb(1:fix(length(perfilb)/2),:);
intra = perfilb(fix(length(perfilb)/2)+1:length(perfilb),:);
 
figure()
hold on
plot (x,y)

as1 = rand(11,1)/1000;  %%
as2 = -rand(11,1)/1000; %% 

perfil = update(perfilb,as1,as2);
x=perfil(:,1);
y=perfil(:,2);
plot (x,y)

%res=analise(perfil)
[ang,cl,cd] = analise(perfil,0);


for i = 1:1
    
i
    
updte=zeros(1,length(as1));
for k = 1:length(as1)
    up1 = 0;
    up2 = 0;
    up3 = 0;
    as1_ = zeros(1,length(as1),1);
    as2_ = zeros(1,length(as1));
    step = as1(k)*1.01;
    as1_(k) = step;
    perfil_ = update(perfilb,as1_+as1,as2);
    [ang_,cl_,cd_] = analise(perfil_,k);
    deriv = (cl_-cl)/(step);
    for a = 1:length(ang_)
        a = ang_(a);
        if a == ang(1)
            up1 = 2*(0.1 - cl_(1))*deriv(1)*0.1;
        elseif a == ang(2)
            up2 = 2*(0.5 - cl_(2))*deriv(2)*0.1;
        elseif a == ang(3)
            up3 = 2*(1 - cl_(3))*deriv(3)*0.1;
        end
    
    end
    
    updte(k) = up2;
    
end
as1_temp = as1+updte;

updte=zeros(1,length(as2));
for k = 1:length(as2)
    up1 = 0;
    up2 = 0;
    up3 = 0;
    as1_ = zeros(1,length(as2),1);
    as2_ = zeros(1,length(as2));
    step = as2(k)*1.1;
    as2_(k) = step;
    perfil_ = update(perfilb,as1,as2_+as2);
    [ang_,cl_,cd_] = analise(perfil_,k);
    deriv = (cl_-cl)/(step);
    for a = 1:length(ang_)
        a = ang_(a);
        if a == ang(1)
            up1 = 2*(0.1 - cl_(1))*deriv(1)*0.1;
        elseif a == ang(2)
            up2 = 2*(0.5 - cl_(2))*deriv(2)*0.1;
        elseif a == ang(3)
            up3 = 2*(1 - cl_(3))*deriv(3)*0.1;
        end
    
    end
    
    updte(k) = up2;
    
end

as1 = as1_temp;
as2 = as2+updte;

perfil = update(perfilb,as1,as2);
x=perfil(:,1);
y=perfil(:,2);
plot (x,y)

% [ang,cl,cd] = analise(perfil,0);

end

function perfil_ = update(perfil,as1,as2)
%%calcula o novo perfil
x=perfil(:,1);
y=perfil(:,2);
extra = perfil(1:fix(length(perfil)/2),:);
intra = perfil(fix(length(perfil)/2)+1:length(perfil),:);
% a's iniciais
t1s = linspace(0.1,0.9999999,11);
t2 = 5;
%as1 = rand(11,1)/1000;
%as2 = -rand(11,1)/1000;
x_ = flip(extra(:,1));
%figure
for cnt = 1:length(t1s)
    t1 = t1s(cnt);
    a1 = as1(cnt);
    %x_ = flip(extra(:,1));
    %y_ = x_;
    y_ = a1 * sin(pi*x_.^(log(0.5)/log(t1))).^t2;
    %plot(x_,y_)
    %hold on
    extra(:,2) = extra(:,2) + y_;
end
x_ = flip(intra(:,1));
for cnt = 1:length(t1s)
    t1 = t1s(cnt);
    a2 = as2(cnt);
    %x_ = flip(extra(:,1));
    %y_ = x_;
    y_ = a2 * sin(pi*x_.^(log(0.5)/log(t1))).^t2;
    %plot(x_,y_)
    %hold on
    intra(:,2) = intra(:,2) + y_;
end
perfil_ = [extra ; intra];
x=perfil_(:,1);
y=perfil_(:,2);
%plot (x,y)
end

%function resultados = analise(perfil)
function [ang, cl, cd] = analise(perfil,k)
    %%faz as analises e diz resultados
    
    %inputs
    Re = 10e5;
    angles = [1,3,5];
    sname = int2str(k) + "save.txt";% nameDat(1) +
    delete(sname)
    %sname = "save.txt";
    iname = "_input.txt";%nameDat(1) + 

    fname = 'analise.txt';%ficheiro com o novo perfil $$$$$$$$$$$$$$$
    pointsX = perfil(:,1);
    pointsY = perfil(:,2);
    
    %% write coordinates
    fid = fopen(fname,'w');
    %fprintf(fid,"%s AIRFOIL \n", "Airfoil");

    % write down all the points
    i = 1;
    while i <= length(pointsX)
        fprintf(fid,"  %f  %f\n",pointsX(i), pointsY(i));
        i = i + 1;
    end
    
    % Close file
    fclose(fid);
    
    fid = fopen(iname,'w');
    
    fprintf(fid,"PLOP\n"); %abrir op????es dos gr??ficos
    fprintf(fid,"G\n\n"); %impedir a sua cria????o
    
    %% Loading the airfoil
    fprintf(fid,"load %s\n",fname);
    fprintf(fid,'\n\n');
    fprintf(fid,'PPAR\n');
    fprintf(fid,"N %i \n", 200);
    fprintf(fid,'\n\n');
    
    %% Find the Cp vs. X plot
    fprintf(fid,'OPER\n');
    fprintf(fid,'Visc\n');
    fprintf(fid,'%i\n',Re);  %reynolds number
    fprintf(fid,'ITER 200\n');
    fprintf(fid,'PACC\n');
    fprintf(fid,"%s \n\n", sname);  %onde se escreve os resultados
    
%     for i=1:(length(config)-config(22))/2  %number of angles
%         fprintf(fid,"Alfa %1.1f \n", config(config(22)-1+2*i));
%     end

    for i=1:length(angles)  %number of angles
        fprintf(fid,"Alfa %1.1f \n", angles(i));
    end
    
%     if config(7+config(21))~=0
%         fprintf(fid,"Alfa %1.1f \n", config(6+config(21)));
%     end
    
    fprintf(fid,'PACC\n');
    fprintf(fid,'\n\n');
    
    % Close file
    fclose(fid);
    
    
    % Initiate the evaluation
    command = 'xfoil.exe < '+ iname + ' &';
    system(command);
    
    disp("pausing")
    pause(3);
    
    %Close all Xfoil unfinished tasks
    system('TASKKILL /F /IM cmd.exe /T ');
    system('TASKKILL /IM xfoil.exe');
    
    
    %%Resultados
    fidsave = fopen(sname,'r');
    dataBuffer = textscan(fidsave,'   %f   %f   %f   %f  %f   %f   %f','HeaderLines',12);
    fclose(fidsave);
    resultados=dataBuffer(:,1:3);
    ang = cell2mat(resultados(1));
    cl = cell2mat(resultados(2));
    cd = cell2mat(resultados(3));
    
    

end


