step = 0.01;
Step = [0.2 , 0.15 , 0.1 , 0.05 , 0.01];
s = 1;
step = Step(s);
% perfil inicial
perfilb=importdata('naca0009.txt', ' ');
perfilb=importdata('NACA2412.txt', ' ');
%perfilb=importdata('E66.txt', ' ');
x=perfilb(:,1);
y=perfilb(:,2);
extra = perfilb(1:fix(length(perfilb)/2),:);
intra = perfilb(fix(length(perfilb)/2)+1:length(perfilb),:);

fname = "resultados.txt";
fid = fopen(fname,'w');

 
figure()
hold on
axis([0 1 -0.15 0.15])
plot (x,y)

ang = [];
n = 0;
d = 0;
%dorso = "teste";
while length(ang)~=3


as1 = rand(20,1)/10 - 0.05;  %%
as1_0 = as1;
as2 = rand(20,1)/10 - 0.05; %%
as2_0 = as2;
%as1 = zeros(40,1)/10;  %%
%as2 = zeros(40,1)/10; %%


perfil = update(perfilb,as1,as2);
x=perfil(:,1);
y=perfil(:,2);
plot (x,y)

%res=analise(perfil)
[ang,cl,cd] = analise(perfil,0,n,d,0);

end

cl_0 = cl;
cd_0 = cd;

cl_id = [0.3652 ; 0.8611 ; 1.1];
%a = 2;


for i = 1:50 %otimizacoes p o mesmo angulo

    
%     if length(cl) == length(cl_id)
%     err = sum((cl_id-cl).^2 + cd);
%     end

    
    %i
    as1_copy = as1;
    as2_copy = as2;


    %[ang,cl,cd] = analise(perfil,0,n,d);

    fprintf(fid, "%d \n", n)
    for j = 1:length(ang)
    fprintf(fid,"%f %f %f \n", ang(j), cl(j), cd(j));
    %fprintf(fid,"%f\n", err);
    end
    fprintf(fid, "------------------\n")



    %loss = (cl(a)-cl_id(a))^2;

    %mkdir (int2str(i))
    for k = 1:length(as1)
        deriv = 0;
        %for a = 1:3
        cd_ = [];
        deriv_ = 0;
        %as1_ = zeros(length(as1),1);
        %as2_ = zeros(length(as2),1);
        %%% as1_ = as1_copy;
        %as1_(k) = as1_(k)*1.1;
        %as1_(max(1,k-3):min(k+3,length(as1_))) = as1_(max(1,k-3):min(k+3,length(as1_)))*1.1;
        %%% as1_(max(1,k-3):min(k+3,length(as1_))) = as1_(max(1,k-3):min(k+3,length(as1_))) + step;
        %as1_(k) = as1_(k) + 0.1e-3;

        %%% perfil_ = update(perfilb,as1_,as2);
        %%% d = 'extra';
        %%% [ang_, cl_, cd_] = analise(perfil_,k,n,d,a);
        
        while length(cd_)<length(cl_id) && s <= length(Step)
           step = Step(s);
           as1_ = as1_copy;
           %as1_(max(1,k-5):min(k+5,length(as1_))) = as1_(max(1,k-5):min(k+5,length(as1_))) + step;
           as1_(k) = as1_(k)+step;
           perfil_ = update(perfilb,as1_,as2);
           d = 'extra';
           [ang_, cl_, cd_] = analise(perfil_,k,n,d,'a');
           s = s+1;
        end
        
        if length(cd)<length(cd_) && length(cd_)==3
           cd = cd_ ;
           perfil = perfil_;
           as1_copy = as1_;
%            as2_copy = as2_;
%             temp = cd;
%             for c = 1:length(cl_id)
%                 if sum(cl_id(c)==cl)>0
%                     cd(c) = temp(cl == cl_id(c));
%                 else
%                     cd(c) = cd_(c);
%                 end
%             end
%             
        end

        for a = 1:2
        if length(ang_) >= a && length(ang) >= a% && ang_(a) == ang(a)
            deriv_ = cd(a)-cd_(a);
            deriv = deriv + deriv_*100*step;
        end
        end
        if length(ang_) >= 3 && length(ang) >= 3 && ang_(3) == ang(3)
            deriv_ = 2 * (1.3 - cl(3)) * (cl_(a)-cl(3));
            deriv = deriv + deriv_*0.01;
        end
        s = 1;
        %deriv = deriv - (mean(cd_) - mean(cd));

        

        as1(k) = as1(k) + deriv;


    end


    for k = 1:length(as2)
        deriv = 0;

        cd_ = [];
        deriv_ = 0;
        %as1_ = zeros(length(as1),1);
        %as2_ = zeros(length(as2),1);
        %%% as2_ = as2_copy;
        %as2_(k) = as2_(k)*1.1;
        %as2_ = as2;
        %as2_(max(1,k-3):min(k+3,length(as2_))) = as2_(max(1,k-3):min(k+3,length(as2_)))*1.1;
        %%% as2_(max(1,k-3):min(k+3,length(as2_))) = as2_(max(1,k-3):min(k+3,length(as2_))) + step;
        %as2_(k) = as2_(k) + 0.05e-3;

        %%% perfil_ = update(perfilb,as1,as2_);
        %%% d = 'intra';
        %%% [ang_, cl_, cd_] = analise(perfil_,k,n,d,a);

        while length(cd_)<length(cl_id) && s <= length(Step)
           step = Step(s);
           as2_ = as2_copy;
           %as2_(max(1,k-5):min(k+5,length(as2_))) = as2_(max(1,k-5):min(k+5,length(as2_))) + step;
           as2_(k) = as2_(k)+step;
           perfil_ = update(perfilb,as1_,as2);
           d = 'intra';
           [ang_, cl_, cd_] = analise(perfil_,k,n,d,a);
           s = s+1;
        end
%         s = 1;
        if length(cd)<length(cd_)
           cd = cd_ ;
           perfil = perfil_;
           as1_copy = as1_;
           as2_copy = as2_;
        end

        for a = 1:2
        if length(ang_) >= a && length(ang) >= a% && ang_(a) == ang(a)
            deriv_ = cd(a)-cd_(a);
            deriv = deriv + deriv_*100*step;
        end
        end
        if length(ang_) >= 3 && length(ang) >= 3 && ang_(3) == ang(3)
            deriv_ = 2 * (1.3 - cl(3)) * (cl_(a)-cl(3));
            deriv = deriv + deriv_*0.01;
        end
        s = 1;
        %deriv = deriv - (mean(cd_) - mean(cd));

        

        as2(k) = as2(k) + deriv;

    end 




    % else
    % loss = mean(cd);
    % 
    % %mkdir (int2str(i))
    % for k = 1:length(as1)
    %     deriv = 0
    %     for a = 1:4%numero de angulos (4º angulo é o Cd)
    %     %as1_ = zeros(length(as1),1);
    %     %as2_ = zeros(length(as2),1);
    %     %as1_(k) = as1(k)*1.2;
    %     %perfil_ = update(perfilb,as1_,as2_);
    %     %[ang_, cl_, cd_] = analise(perfil_,k);
    %     as1_ = as1_copy;
    %     as1_(max(1,k-3):min(k+3,length(as1_))) = as1_(max(1,k-3):min(k+3,length(as1_)))*1.2;
    %     perfil_ = update(perfilb,as1_,as2_copy);
    %     d = 'extra';
    %     [ang_, cl_, cd_] = analise(perfil_,k,n,d);
    %     
    %     if length(cl)<length(cl_)
    %        cl = cl_ ;
    %        cd = cd_;
    %     end
    %     
    %     if length(cd_)>0 
    %     %if length(ang_) == a && ang_(a) == ang(a)
    %     loss_ = mean(cd_);
    %     %as1(k) = as1(k)*(loss<loss_)*0.8 + as1_(k)*(loss>loss_);
    %     deriv_ = 2 * loss * (loss_ - loss);
    %     if isnan(loss_) == 1 || isnan(loss) == 1 || isnan(deriv) == 1
    %     deriv_ = 0;
    %     end
    %     
    %     deriv = deriv + deriv_;
    %     
    %     end
    %     end
    %     
    %     as1(k) = as1(k) + deriv*200;



    % for k = 1:length(as2)
    %     %as1_ = zeros(length(as1),1);
    %     %as2_ = zeros(length(as2),1);
    %     %as2_(k) = as1(k)*1.2;
    %     %perfil_ = update(perfilb,as1_,as2_);
    %     %[ang_, cl_, cd_] = analise(perfil_,k);
    %     as2_ = as2_copy;
    %     as2_(max(1,k-3):min(k+3,length(as2_))) = as2_(max(1,k-3):min(k+3,length(as2_)))*1.2;
    %     perfil_ = update(perfilb,as1_copy,as2_);
    %     d = 'intra';
    %     [ang_, cl_, cd_] = analise(perfil_,k,n,d);
    %     if length(cl)<length(cl_)
    %        cl = cl_ ;
    %        cd = cd_;
    %     end
    %     
    %     if length(cd_)>0 
    %     %if length(ang_) >= a && ang_(a) == ang(a)
    %     loss_ = mean(cd_);
    %     %as2(k) = as2(k)*(loss<loss_)*0.8 + as2_(k)*(loss>loss_);
    %     deriv = 2 * loss * (loss_ - loss);
    %     if isnan(loss_) == 0 && isnan(loss) == 0 && isnan(deriv) == 0
    %     as2(k) = as2(k) + deriv*200;
    %     end
    %     end
    %     
    % end


    perfil = update(perfilb,as1,as2);
    x = perfil(:,1);
    y = perfil(:,2);
    plot (x,y)

    d = 'novo';
    %res=analise(perfil)
    [ang,cl,cd] = analise(perfil,0,n,d, 0);

    n = n + 1;
    %step = 0.1 + 0.1*(length(cl) == 3);

end

%end
%n = n + 1


%fclose(fname);

function perfil_ = update(perfil,as1,as2)
%%calcula o novo perfil

x=perfil(:,1);
y=perfil(:,2);
extra = perfil(1:fix(length(perfil)/2),:);
intra = perfil(fix(length(perfil)/2)+1:length(perfil),:);
% a's iniciais
t1s = linspace(0.1,0.9999999,length(as1));
t2 = 2;
%as1 = rand(11,1)/1000;
%as2 = -rand(11,1)/1000;
%x_ = flip(extra(:,1));
x_ = extra(:,1);
%figure
for cnt = 1:length(t1s)
    t1 = t1s(cnt);
    a1 = as1(cnt);
    %x_ = flip(extra(:,1));
    %y_ = x_;
    y_ = a1 * 0.01* sin(pi*x_.^(log(0.5)/log(t1))).^t2;
    %plot(x_,y_)
    %hold on
    extra(:,2) = extra(:,2) + y_;
end
%x_ = flip(intra(:,1));
x_ = (intra(:,1));
for cnt = 1:length(t1s)
    t1 = t1s(cnt);
    a2 = as2(cnt);
    %x_ = flip(extra(:,1));
    %y_ = x_;
    y_ = a2 * 0.01* sin(pi*x_.^(log(0.5)/log(t1))).^t2;
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
function [ang, cl, cd] = analise(perfil,k,n,d,a)
    %%faz as analises e diz resultados
    g = pwd;
    %inputs
    Re = 3.5e5;
    angles = [1,8,10];
    cls = [0.2665 ; 0.5750];
    %if d == 1
    %    dorso = "extra"
    %elseif d == 2
    %    dorso == "intra"
    %elseif d == 0
        %não fazer nada
    %end    
    
        
    
    curfold = int2str(n) + "pasta"; %current folder - nome da pasta atual
    
    if not(isfolder(curfold))
        disp('ok')
        mkdir(g, [curfold])%cria uma pasta em ../geneticos3 com
        %o nome definido na variável curfold, que depende da iteração do
        %perfil em q o programa se encontra
    end
    

    
    
    %mkdir ../geneticos3 foldname
    
    sname = "a" + a + "_" + d + "_" + int2str(k) + "save.txt";% nameDat(1) +  %%% resultados
    %destpath = (['C:\Users\mbv50\OneDrive - Universidade de Lisboa\Faculdade - IST\AeroTéc - OAT\ACC-2022\Programa Genético\Código\git_win\geneticos3\', num2str(curfold)])
    %movefile(sname,destpath) - movido p o fim da funcao
    destpath = g + "\" + curfold;
    %destpath = ([destpath, num2str(curfold)]);
    
    
    %delete(sname)
    %sname = "save.txt";
    
    iname = "_input.txt";%nameDat(1) + 

    fname =   d + "_" + int2str(k) + "analise.txt";%ficheiro com o novo perfil $$$$$$$$$$$$$$$
    %movefile(fname,destpath) - movido p o fim da funcao
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
    
    fprintf(fid,"PLOP\n"); %abrir opções dos gráficos
    fprintf(fid,"G\n\n"); %impedir a sua criação
    
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
    
    
    for i=1:length(cls)  %number of angles
        %fprintf(fid,"Alfa %1.1f \n", 10);
        fprintf(fid,"Cl %1.4f \n", cls(i));
    end
    fprintf(fid,"Alfa %1.1f \n", 12);
    
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
    pause(3 + (k==0)*2);
    
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
    
    movefile(sname,destpath)%coloca o ficheiro sname na pasta correspondente à iteração do perfil em que o programa está
    movefile(fname,destpath)%coloca o ficheiro fname na pasta correspondente à iteração do perfil em que o programa está

end


