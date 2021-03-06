%{ 
function release2(varargin)
%the code has the option to be called from another one. If the correct
%input is given it will run. if not it will assume to be running in solo
%mode, asking for the settings.

%checking number of inputs auxiliary calculations
setingsNumber = 21; %acconting for the values that are not angles nor weights
modeNumber = nargin - setingsNumber;

%% define options based on operation mode
if rem(modeNumber,2) == 0 && modeNumber > 0
    %use the input configuration
    %
    % varargin{1} = best individual lenght
    % varargin{2} = generation limit
    % varargin{3} = number of offspring 1
    % varargin{4} = Generation limit 1
    % varargin{5} = Max thickness Bonus percentage 1
    % varargin{6} = CL angle to check 1
    % varargin{7} = CL value to ckeck 1
    % varargin{8} = generational change (1 = yes 0 = no)
    % varargin{9} = generation change limit 
    % varargin{10} = number of offspring 2
    % varargin{11} = Generation limit 2
    % varargin{12} = Max thickness Bonus percentage 2
    % varargin{13} = CL angle to check 2
    % varargin{14} = CL value to ckeck 2
    % varargin{15} = Number of Points
    % varargin{16} = Re Number
    % varargin{17} = Type of Point Deduction
    % varargin{18} = Enhanced Aproximation
    % varargin{19} = Starting Time
    % varargin{20} = Maximum allowed time wait for profiles
    % varargin{21} = phase
    
    % varargin{22} = angle 1
    % varargin{23} = angle weight   
    % etc
    
    mode = "input";
    
    config = zeros(1,nargin);
    
    for n = 1:nargin
        config(1,n) = varargin{n};
    end
    
elseif nargin == 0
    %} 

setingsNumber = 24; %acconting for the values that are not angles nor weights ------------
%modeNumber = nargin - setingsNumber; %-------------------------------

    %% Configuration UI
    %create a new configuration based on ui
    
    %UI's setup 
    fig = uifigure('Name','Program Inputs','NumberTitle','off', 'Position', [50 50 800 480]); %criar janela
    
    
    %Run button
    btn = uibutton(fig,'push','Text','Run',...
        'Position',[450, 30, 100, 30],...
        'ButtonPushedFcn', @(btn,event) uiresume(fig)); 
    
    labelRun1 = uilabel(fig,'Text','Fields marked with "*" are','FontSize',12,...
        'FontColor','k','Position',[600, 35, 400, 40]);
    
    labelRun2 = uilabel(fig,'Text','removed with a "0" value','FontSize',12,...
        'FontColor','k','Position',[600, 20, 400, 40]);
    
    
    
    %Stop criteria
    labelStop = uilabel(fig,'Text','Stopping Criteria','FontSize',20,...
        'FontWeight','bold','FontColor','k','Position',[100, 440, 400, 40]);
    
    %Best individual
    labelGenBestMax = uilabel(fig,'Text','Best individual limit*','FontSize',12,...
        'FontColor','k','Position',[100, 420, 400, 40]);
    
    genBestMax = uieditfield(fig,'numeric','Position',[100 405 100 25], ...
        'RoundFractionalValues','on','Value', 50);
    
    genBestMax.Limits = [0 1000];    
    
    %Max Generation
    labelGenMax = uilabel(fig,'Text','Generation Limit*','FontSize',12,...
        'FontColor','k','Position',[100, 375, 400, 40]);
    
    genMax = uieditfield(fig,'numeric','Position',[100 360 100 25], ...
        'RoundFractionalValues','on','Value', 400);
    
    genMax.Limits = [1 10000];
    
    
    
    %Specific Criteria
    labelSpecific = uilabel(fig,'Text','Specific Criteria','FontSize',20,...
        'FontWeight','bold','FontColor','k','Position',[275, 440, 400, 40]);
    
    %Number of points
    labelNPoints= uilabel(fig,'Text','Number of Points','FontSize',12,...
        'FontColor','k','Position',[275, 420, 400, 40]);
    
    NPoints = uieditfield(fig,'numeric','Position',[275 405 100 25], ...
        'RoundFractionalValues','on','Value', 11);
    
    NPoints.Limits = [4 100];    
    
    %Re Number
    labelReNumber = uilabel(fig,'Text','Re Number','FontSize',12,...
        'FontColor','k','Position',[275, 375, 400, 40]);
    
    ReNumber = uieditfield(fig,'numeric','Position',[275 360 100 25], ...
        'RoundFractionalValues','on','Value', 400000);
    
    ReNumber.Limits = [100 1000000000000];
    
    %Type of Point Deduction
    %labelTypeDeduction = uilabel(fig,'Text','Type of Point Deduction','FontSize',12,...
    %    'FontColor','k','Position',[275, 330, 400, 40]);
    
    typeDeduction = uiswitch(fig,'rocker','Items',{'Custom','Xfoil'},'Position',[330 280 30 60]);
    
    
    %Enhanced Aproximation
    %labelEnhanced = uilabel(fig,'Text','Enhanced Aproximation','FontSize',12,...
    %    'FontColor','k','Position',[275, 285, 400, 40]);
    
    enhanced = uiswitch(fig,'toggle','Items',{'No','Yes'},'Position',[280 280 30 60]);
    
    %Maximum calculation time
    labelMaxTime = uilabel(fig,'Text','Maximum generation time(s)','FontSize',12,...
        'FontColor','k','Position',[275, 235, 400, 40]);
    
    maxTime = uieditfield(fig,'numeric','Position',[275 220 100 25], ...
        'RoundFractionalValues','off','Value', 12);
    
    maxTime.Limits = [0 60];
    
    
    
    %Options
    labelStop = uilabel(fig,'Text','Options','FontSize',20,...
        'FontWeight','bold','FontColor','k','Position',[100, 300, 400, 40]);
    
    %Number of offspring
    labelNOfspring = uilabel(fig,'Text','Number of Offspring','FontSize',12,...
        'FontColor','k','Position',[100, 280, 400, 40]);
    
    NOfspring = uieditfield(fig,'numeric','Position',[100 265 100 25], ...
        'RoundFractionalValues','on','Value', 20);
    
    NOfspring.Limits = [0 1000];    
    
    %Father racio
    labelFatherRacio = uilabel(fig,'Text','Paternal Racio','FontSize',12,...
        'FontColor','k','Position',[100, 235, 400, 40]);
    
    FatherRacio = uieditfield(fig,'numeric','Position',[100 220 100 25], ...
        'RoundFractionalValues','off','Value', 0.50);
    
    FatherRacio.Limits = [0 1];

    %Max thickness bonus
    labelThickBonus = uilabel(fig,'Text','Max Thickness Bonus (percentage)*','FontSize',12,...
        'FontColor','k','Position',[100, 190, 400, 40]);
    
    thickBonus = uieditfield(fig,'numeric','Position',[100 175 100 25], ...
        'RoundFractionalValues','off','Value', 9.5);
    
    thickBonus.Limits = [0 100];
    
    %CL Check angle
    labelThickBonusAlfa = uilabel(fig,'Text','CL angle to check','FontSize',12,...
        'FontColor','k','Position',[100, 145, 400, 40]);
    
    CLAlfa = uieditfield(fig,'numeric','Position',[100 130 100 25], ...
        'RoundFractionalValues','off','Value', 2.5);
    
    CLAlfa.Limits = [-10 20];
    
    %CL Value Check
    labelCLBonus = uilabel(fig,'Text','CL value to check*','FontSize',12,...
        'FontColor','k','Position',[100, 100, 400, 40]);
    
    CLBonus = uieditfield(fig,'numeric','Position',[100 85 100 25], ...
        'RoundFractionalValues','off','Value', 0);
    
    CLBonus.Limits = [0 5];
    
    %Evaluation
    labelEvalAng = uilabel(fig,'Text','Angles and its percentage','FontSize',12,...
        'FontColor','k','Position',[100, 55, 400, 40]);
    
    evalAng = uieditfield(fig,'Position',[100 40 300 25], ...
        'Value', "0,0.3;2,0.2;4,0.2;5,0.3");
    
    %stall detection
    %stall angle to check
    labelStalAngle = uilabel(fig,'Text','stall angle to check','FontSize',12,...
        'FontColor','k','Position',[300, 145, 400, 40]);
    
    StalAngle = uieditfield(fig,'numeric','Position',[300 130 100 25], ...
        'RoundFractionalValues','off','Value', 12);
    
    StalAngle.Limits = [-30 30];
    
    %Allowed Stal variance
    labelAllowedStalVariance = uilabel(fig,'Text','Allowed Variance*','FontSize',12,...
        'FontColor','k','Position',[300, 100, 400, 40]);
    
    AllowedStalVariance = uieditfield(fig,'numeric','Position',[300 85 100 25], ...
        'RoundFractionalValues','off','Value', 0.2);
    
    AllowedStalVariance.Limits = [0 1];
    
    %Generational change
    labelGenChange = uilabel(fig,'Text','Generational Change','FontSize',20,...
        'FontWeight','bold','FontColor','k','Position',[450, 440, 400, 40]);
    
    genChange = uiswitch(fig,'Items',{'No','Yes'},'Position',[470 410 60 90]);

    %Generational change position
    labelGenChangePos = uilabel(fig,'Text','Generation Change Limit','FontSize',12,...
        'FontColor','k','Position',[450, 375, 400, 40]);
    
    genChangePos = uieditfield(fig,'numeric','Position',[450 360 100 25], ...
        'RoundFractionalValues','on','Value', 150);  
    
    genChangePos.Limits = [0 10000];
    
    
    
    %Options 2
    labelStop2 = uilabel(fig,'Text','Options','FontSize',20,...
        'FontWeight','bold','FontColor','k','Position',[450, 300, 400, 40]);
    
    %Number of offspring 2
    labelNOfspring2 = uilabel(fig,'Text','Number of Offspring','FontSize',12,...
        'FontColor','k','Position',[450, 280, 400, 40]);
    
    NOfspring2 = uieditfield(fig,'numeric','Position',[450 265 100 25], ...
        'RoundFractionalValues','on','Value', 20);
    
    NOfspring2.Limits = [0 1000];    
    
    %Father racio 2
    labelFatherRacio2 = uilabel(fig,'Text','Generation Limit','FontSize',12,...
        'FontColor','k','Position',[450, 235, 400, 40]);
    
    FatherRacio2 = uieditfield(fig,'numeric','Position',[450 220 100 25], ...
        'RoundFractionalValues','off','Value', 0.50);
    
    FatherRacio2.Limits = [0 1];

    %Max thickness bonus 2
    labelThickBonus2 = uilabel(fig,'Text','Max Thickness Bonus (percentage)','FontSize',12,...
        'FontColor','k','Position',[450, 190, 400, 40]);
    
    thickBonus2 = uieditfield(fig,'numeric','Position',[450 175 100 25], ...
        'RoundFractionalValues','off','Value', 9.5);
    
    thickBonus2.Limits = [0 100];
    
    %CL Check angle 2
    labelThickCLAlfa2 = uilabel(fig,'Text','CL angle to check','FontSize',12,...
        'FontColor','k','Position',[450, 145, 400, 40]);
    
    CLAlfa2 = uieditfield(fig,'numeric','Position',[450 130 100 25], ...
        'RoundFractionalValues','off','Value', 1.5);
    
    CLAlfa2.Limits = [-10 20];
    
    %CL Value Check 2
    labelCLBonus2 = uilabel(fig,'Text','CL value to check','FontSize',12,...
        'FontColor','k','Position',[450, 100, 400, 40]);
    
    CLBonus2 = uieditfield(fig,'numeric','Position',[450 85 100 25], ...
        'RoundFractionalValues','off','Value', 1.1);
    
    CLBonus2.Limits = [0 5];
    
    
    
    %% Value interpretation
    uiwait(fig); 
    
    %Check Generation Change Switch
    if strcmp(genChange.Value,'Yes')
        genChangeValue = 1;
    else
        genChangeValue = 0;
    end
    
    %Check Generation Change Switch
    if strcmp(typeDeduction.Value,'Xfoil')
        typeDeductionValue = 1;
    else
        typeDeductionValue = 0;
    end
    
    %Check Generation Change Switch
    if strcmp(enhanced.Value,'Yes')
        enhancedValue = 1;
    else
        enhancedValue = 0;
    end
    
    evalAngSplit = split(evalAng.Value,";");
    config = zeros(1,setingsNumber+2*length(evalAngSplit));

    %evalAngSplit2
    config(1,1) = genBestMax.Value;
    config(1,2) = genMax.Value;
    config(1,3) = NOfspring.Value;
    config(1,4) = FatherRacio.Value;
    config(1,5) = thickBonus.Value;
    config(1,6) = CLAlfa.Value;
    config(1,7) = CLBonus.Value;
    config(1,8) = genChangeValue;
    config(1,9) = genChangePos.Value;
    config(1,10) = NOfspring2.Value;
    config(1,11) = FatherRacio2.Value;
    config(1,12) = thickBonus2.Value;
    config(1,13) = CLAlfa2.Value;
    config(1,14) = CLBonus2.Value;    
    config(1,15) = NPoints.Value;
    config(1,16) = ReNumber.Value;
    config(1,17) = typeDeductionValue;
    config(1,18) = enhancedValue;
    config(1,19) = str2double(datestr(now, 'ddmmyyHHMM')); %starting time
    config(1,20) = maxTime.Value;
    config(1,21) = 0; %phase
    config(1,22) = setingsNumber;
    config(1,23) = AllowedStalVariance.Value;
    config(1,24) = StalAngle.Value;
    
    percentageCheck = 0;
    for n = 1:length(evalAngSplit)
        evalAngSplit2 = split(evalAngSplit{n},",");
        config(1,setingsNumber-1+2*n) = str2double(evalAngSplit2{1});
        config(1,setingsNumber+2*n) = str2double(evalAngSplit2{2});
        percentageCheck = percentageCheck + str2double(evalAngSplit2{2});
        
    end
    
    if percentageCheck ~= 1
        disp("warning percentage weights do not amount to 1")
    end
    close(fig);

    %{
else
    %if the number of input arguments is wrong
    disp("incorrect number of input arguments");
    return 
    
end
%}
%% input gathering and population starting
population = inputProfiles(config);

%% Progres Bar
%{
f = waitbar(0,'Estimated Progress','Name','Maias but in English',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);
%}

%% Main loop
gen = 0;

while true
    %Loop Regulators
    gen = gen + 1
    if gen == config(9)
        config(21) = 7;
    end
    
    %criar configs e ficheiros de perfil
    population = createinput(population, config);
    
    %xfoil paralelo
    population = parallel(population, config);
    
    %rank
    [population,best] = selec(population,config);
    
    %Save Progress
    saveData(population,config);
    
    %Stopping Check
    if (gen > config(2) && config(2)>0) || (population(best).gen < gen-config(1) && config(1)>0)
        break
    end
    
    %{
    if getappdata(f,'canceling')
        break
    end
    
    % Update waitbar and message
    if config(2)>0
        waitbar(gen/config(2),f,sprintf('Estimated Progress'))
    end
    %}
    
    %Generating offsprings
    population = recomb(population,config,gen);

end

%close main ui
%delete(f);

%Save the variables for future use
save('popu.mat')

%end %end complete function 



%% Auxiliary Functions
function population = inputProfiles(config)
    %function that iniciates the population and inserts the initial
    %profiles, converting them into the reduced points form
    
    
    listing = dir("initialAirfoils");
    population = struct;
    ignored = 0;
    
    for n = 1:length(listing)
        
        %check if the file type is adequate
        fileCheck = split(listing(n).name,".");
        
        if strcmp(fileCheck(length(fileCheck)),'txt') || strcmp(fileCheck(length(fileCheck)),'dat')
            %setup the adequate candidate
            population(n-ignored).name = join(['createdAirfoils/' listing(n).name],'');
            population(n-ignored).gen = 0;
            population(n-ignored).n = n-ignored;
            
            %colecting airfoil data
            local = join(["initialAirfoils/" listing(n).name],"");
            data=importdata(local(1), ' ', 1);
            x=data.data(:,1);
            y=data.data(:,2); 
            %p=length(x); %number of points in the airfoil file
            
            
            %separate the airfoil
            [~,j]=min(x); %x from 0 to 1 -> x=0 ?? o BA
            x_extra=x(1:j-1);
            x_intra=x(j:end); %inclui BA em x=0

            nfi=config(1,15); %Number of control points
            fi=zeros(1,nfi);
            xfi=zeros(1,nfi);
            yfi=zeros(1,nfi);

            for i=1:nfi

                fi(i)=(i-1)*360/(nfi-1); % nfi pontos de controlo espa??ados entre 0 e 360
                xfi(i)=0.5*(1+cos( pi*(fi(i))/180 )); %coordenada x correspondente

                if fi(i)<=180 %intradorso
                    dif=x_intra-xfi(i);
                    [~,indice]=min(abs(dif));
                    yfi(i)=y(indice+j-1);
                    %xfi(i)=x(indice+j-1); %ajuste do x
                else %extradorso
                    dif=x_extra-xfi(i);
                    [~,indice]=min(abs(dif));
                    yfi(i)=y(indice);
                    %xfi(i)=x(indice);  %ajuste do x
                end


            end

            population(n-ignored).control(:,1)=xfi;
            population(n-ignored).control(:,2)=yfi;

        else
            %if the found document is not adequate the it is skiped
            ignored = ignored + 1;
            
        end
    end
end


function population = createinput(population, config)
if population(length(population)).gen==0
    base = 1;
else
    base = length(population)-config(3+config(21));
end

for n = base:length(population)
    nameDat = split(population(n).name,".");
    fname = population(n).name;
    sname = nameDat(1) + "save.txt";
    iname = nameDat(1) + "_input.txt";
    
    %% Coordinates
    %Fun????o adaptada do DrawFoil que recebe o vetor de pontos de controlo e
    %obt??m as coordenadas x,y do perfil. Esta ser?? a fun????o que estar?? a
    %funcionar no loop do algoritmo gen??tico. recebe o 'DNA' do perfil, que
    %consiste num vetor em que cada entrada corresponde a um angulo fi e tem um
    %valor de y/c correspondente. Por exemplo: a entrada numero 1 do vetor
    %corresponde a fi=360, a entrada 2 corresponde a (2-1)*360/(nfi-1), com
    %nfi= legnth(v), etc. E ?? entrada 1 corresponde a coordenada y de fi=360??
    %
    %O perfil inicial tem de ser definido com corda horizontal e em y=0;
    %bordo de ataque para x=0 e bordo de fuga em x=1
    %
    %os pontos de controlo obtidos est??o em fun????o de um ??ngulo fi e igualemnte
    %espa??ados. o ??ngulo fi ?? medido a partir de um referencial centrado no
    %meio da corda e ?? medido em rela????o ao BF no sentido hor??rio
    %
    %no fim as coordenadas do 'novo' perfil s??o obtidas atrav??s de um spline
    %c??bico que passa pelos pontos do gr??fico de y em fun????o de fi
    %HF
    
    %Nfi ?? o numero de pontos de controlo; para garantir que h?? sempre um ponto
    %de controlo no BA (x,y)=(0,0), ?? necess??rio colocar um numero impar
    v = population(n).control;
    X1 = v(:,1);
    X2 = v(:,2);
    v = X2;
    nfi=length(v);
       
    %check for aproximation type
    if config(17) == 0
        p=100; %use 200 points to describe the profile
    else
        p=config(15); %use the xfoil aproximation
    end
    
    %{ 
    %check for enhanced aproximation
    if config(18) == 1
        %enhanced aproximation
    else
        %regular aproximation
    end
    %}
    
    %alocate memory
    fi=zeros(1,nfi);
    xfi=zeros(1,nfi);
    yfi=zeros(1,nfi);
    fic=zeros(1,p);
    yf=zeros(1,p);
    x=zeros(1,p);
    
    %aloca????o de y v2
    xfi = X1;
    yfi = X2;
    
    
    %parte do henrique para calcular y, que ?? agora, em principio, desnecess??ria
    for i=1:nfi
        
        fi(i)=(i-1)*360/(nfi-1); % nfi pontos de controlo espa??ados entre 0 e 360
        xfi(i)=0.5*(1+cos( pi*(fi(i))/180 )); %coordenada x correspondente
        %yfi(i)=v(i);
        
    end
    
    
    xx=linspace(1,360,p);
    %yy=spline(fi,yfi,xx); %spline c??bico para obter coordenadas y do perfil
    %yy=pchip(fi,yfi,xx); %spline c??bico para obter coordenadas y do perfil
    yy=makima(fi,yfi,xx); %spline c??bico para obter coordenadas y do perfil
    
    %obter coordenada x
    for i=1:p
        
        fic(i)=(i-1)*360/(p-1); % nfi pontos de controlo espa??ados entre 0 e 360
        pointsX(i)=0.5*(1+cos( pi*(fic(i))/180 )); %coordenada x correspondente
        yf(i)=yy(i);
        
    end
    
    pointsY=flip(yf);
    
    %plot(pointsX,pointsY)
    
    axis([0 1 -0.2 0.2])
    hold on
    %pointsY = pointsY./(1-pointsX*0.99);
    pointsY = smoothdata(pointsY,'gaussian');
%     points2Y = smoothdata(pointsY,'gaussian');
%     for i=1:length(pointsY)
%         %if pointsX(i) < 0.5 % || i>=length(pointsY)*3/4 
%         pointsY(i) = (pointsX(i)) * pointsY(i) + (1-pointsX(i))*points2Y(i); 
%         %else
%          %   pointsY(i) = (1-pointsX(i)) * pointsY(i) + (pointsX(i))*points2Y(i);
%             %pointsY(i) = (1-pointsX(i)) * pointsY(i) + (pointsX(i))*points2Y(i);
%         %end
%         
%     end    
    
    %pointsY = min(1,-10*pointsX+10).*pointsY;
    
    error = (pointsY(1)-pointsY(length(pointsY)));
%     pointsY(length(pointsY)) = 0;
%     pointsY(1) = 0;
%     pointsY = (min(1,-2*pointsX+2).*pointsY);%*0.5 + pointsY*0.5;
    %pointsY = pointsY.*(1-(3.^pointsX-1)/(3^1-1));
    %yy(1) = yy(1)-error;
    %yy(200) = yy(200)+error;


    %fechar
%     for i=1:50
%         pointsY(i) = pointsY(i)-error/(1+(i-1)/10)*0.5;
%         pointsY(length(pointsY)+1-i) = pointsY(length(pointsY)+1-i)+error/(1+(i-1)/10)*0.5;
%     end
    
%     
%     pointsY = pointsY-((pointsY(50)+pointsY(51))*0.5).*(1-pointsX);
%     pointsY = pointsY-(pointsY(1).*(pointsX));
%     
%     for i=1:50
%         pointsY(i) = pointsY(i)-error/(1+(i-1)/10)*0.5;
%         pointsY(length(pointsY)+1-i) = pointsY(length(pointsY)+1-i)+error/(1+(i-1)/10)*0.5;
%     end


%     erro_c = pointsY(1);
%     erro_b = pointsY(length(pointsY));
%     cima = pointsY(1:length(pointsY)/2);
%     baixo = pointsY(length(pointsY)/2+1:length(pointsY));
%     cima = cima - erro_c*(1-pointsX(1:length(pointsY)/2));
%     baixo = baixo - erro_b*(1-pointsX(length(pointsY)/2+1:length(pointsY)));
%     %pointsX = [pointsX(1:length(pointsY)/2) , 0, pointsX(length(pointsY)/2+1:length(pointsY))];
%     pointsY = [cima, baixo];
%     %thick = max(abs(cima-baixo));
%     %population(n).thick = thick;
%     %if thick <= 0.095
%      %  pointsY = pointsY * (0.10/thick); 
%     %end
    
    plot(pointsX,pointsY)
    
    %no ciclo 'for' procuram-se as coordenadas x e y correspondentes ao fi do
    %ponto de controlo. a coordenada x ?? obtida atr??ves de fi, diretamente.
    %Depois procura-se a entrada do vetor de coordenadas correspondente ao x
    %obtido (xfi). Para isso procura-se o m??nimo do m??dulo do vetor da
    %diferen??a entre as coordenadas x e xfi. o indice do minimo corresponde ??
    %posi????o do y correspondente no vetor de coordenadas y.
    %aten????o que ?? necess??rio procurar os pontos correspondentes para o
    %intradorso se fi<= 180?? e no extradorso se fi>180??
    
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
    
    %% Delete files if they exist
    if (exist(iname,'file'))
        delete(iname);
    end
    if (exist(sname,'file'))
        delete(sname);
    end
    
    
    %%write file intructions
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
    fprintf(fid,'%i\n',config(16));
    fprintf(fid,'ITER 200\n');
    fprintf(fid,'PACC\n');
    fprintf(fid,"%s \n\n", sname);
    for i=1:(length(config)-config(22))/2
        fprintf(fid,"Alfa %1.1f \n", config(config(22)-1+2*i));
    end
    
    if config(7+config(21))~=0
        fprintf(fid,"Alfa %1.1f \n", config(6+config(21)));
    end
    
    fprintf(fid,'PACC\n');
    fprintf(fid,'\n\n');
    
    % Close file
    fclose(fid);
    
    % Initiate the evaluation
    command = 'xfoil.exe < '+ iname + ' &';
    system(command);
    
end
    disp("pausing")
    pause(config(20));
    
    %Close all Xfoil unfinished tasks
    
    system('TASKKILL /F /IM cmd.exe /T ');
    system('TASKKILL /IM xfoil.exe');
end


function population = parallel(population, config)
if population(length(population)).gen==0
    base = 1;
else
    base = length(population)-config(3+config(21));
end



for n = base:length(population)
    nameDat = split(population(n).name,".");
    sname = nameDat(1) + "save.txt";
    f = dir(sname) ;

    if isfile(sname) && f.bytes ~= 0
        fidsave = fopen(sname,'r');
        dataBuffer = textscan(fidsave,'   %f   %f   %f   %f  %f   %f   %f','HeaderLines',12);
        fclose(fidsave);
        delete(sname)
        
        if config(7+config(21)) == 0
            NAngles = (length(config)-config(22))/2;
        else
            NAngles = ((length(config)-config(22))/2)+1;
        end
        
        if length(dataBuffer{1})>= NAngles
            chord(1)=population(n).control(2,2)-population(n).control(10,2);
            chord(2)=population(n).control(3,2)-population(n).control(9,2);
            chord(3)=population(n).control(4,2)-population(n).control(8,2);
            chord(4)=population(n).control(5,2)-population(n).control(7,2);
            chord_max=max(abs(chord));
            
            CLCD=0;
            alfa = zeros ((length(config)-config(22))/2,1);
            CL = zeros ((length(config)-config(22))/2,1);
            CD = zeros ((length(config)-config(22))/2,1);
            CLCD_aux = zeros ((length(config)-config(22))/2,1);
            
            for j=1:(length(config)-config(22))/2
                % Separate Cp data
                alfa(j)  = dataBuffer{1}(j);
                CL(j) = dataBuffer{2}(j);
                CD(j) = dataBuffer{3}(j);
                CLCD_aux(j) = CL(j)/CD(j);
                CLCD = CLCD + CLCD_aux(j)*config(1,2*j+config(22));
            end
            
            if config(5+config(21)) ~= 0 && (chord_max < (config(5+config(21))/100))
                CLCD=CLCD*(abs((chord_max/((config(5+config(21))/100)^2-0.1))));
            end
            
            if config(7+config(21))~=0 && (dataBuffer{2}(j+1)<config(7+config(21)))
                CLCD=0;
            end
            
            if CLCD > 180
                CLCD = 0;
            end
            

            
            sig_cl1 = 1/(1+exp(2-15*CL(1)));
            sig_cl2 = 1/(1+exp(6-14*CL(2)));
            sig_cl3 = 1/(1+exp(9.9785-10.7527*CL(3)));
            population(n).fitness = sig_cl1 * sig_cl2 * sig_cl3;
            % 2,0.3;1.7,0.3;9.6,0.4
            
            %population(n).fitness = CLCD;
            %population(n).fitness = (CLCD * 0.1) * sum(abs(population(n).control(:,2)))*(population(n).control(18)*15-population(n).control(16)*8+population(n).control(13));
            population(n).CL = CL;
            population(n).CD = CD;
            %population(n).CL = espessura;
            %population(n).CD = espessura_c;
            population(n).CLCD = CLCD_aux;
        else
            delete(sname)
            % Null data
            population(n).fitness = 0;
            population(n).CL = 0;
            population(n).CD = 0;
            population(n).CLCD = 0;
        end
    else
        delete(sname)
        % Null data
        population(n).fitness = 0;
        population(n).CL = 0;
        population(n).CD = 0;
        population(n).CLCD = 0;
    end

end
end


function saveData(population, config)
    day = floor(config(19)/10^8);
    month = floor((config(19)/10^6)-100*day);
    year = floor((config(19)/10^4)-10000*day-100*month);
    hour = floor((config(19)/10^2)-1000000*day-10000*month-100*year);
    minutes = floor((config(19))-100000000*day-1000000*month-10000*year-100*hour);
    nameOfFileA = join(["gen_last-" day "_" month "_" year "__" hour "," minutes '.csv'],"");
    nameOfFileB = join(["gen_history-" day "_" month "_" year "__" hour "," minutes '.csv'],"");
    %fid = fopen(nameOfFile,"w");
    %close(fid);
    if length(population)>1
        writetable(struct2table(population),nameOfFileA)
        %writetable
    end
    
end


function [population,best_gen] = selec(population,config)

num_individuos=length(population);
rank_excluir_zeros=[];

%Ver se tem zeros e excluir
for i=1:num_individuos
    if population(i).fitness==0
        rank_excluir_zeros = [rank_excluir_zeros i];
    end
end

for k = 1:length(rank_excluir_zeros)
    nameDat = split(population(rank_excluir_zeros(k)).name,".");
    iname = nameDat(1) + "_input.txt";
    delete(population(rank_excluir_zeros(k)).name);
    delete(iname);
end

population(rank_excluir_zeros)=[];

%Organizar indiv??duos segundo fitness
a = rand();
if a<= 0.8
    b=1;
elseif a>= 0.9
    b=2;
else
    b=3;
end
b=1;
num_individuos=length(population);
rank=zeros(num_individuos,1);
fitness_1=zeros(num_individuos,1);
fitness=zeros(num_individuos,1);
aux=0;
soma=0;


for i=1:num_individuos
    fitness(i)=population(i).fitness(b);
    if fitness(i) <= 0
        fitness(i) = 0;

    end
end

for i=1:num_individuos
    for j=1:num_individuos
        if fitness(j)>=aux
            aux=fitness(j);
            index=j;
        end
    end

    fitness_1(i)=fitness(index);
    fitness(index)=-1;
    %populacao(index).fitness=0;
    rank(i)=index;
    aux=0;
end
best_gen=rank(1);

%Sele????o--------------
if population(length(population)).gen<=3
    num_individuos_select = round(num_individuos*0.8);
else
    num_individuos_select=round(config(4+config(21))*num_individuos);
end
non_select=num_individuos-num_individuos_select;
rank_excluir=zeros(non_select,1);

%Rank relativo----------------
for i=1:num_individuos_select
soma=soma+fitness_1(i);
end

for i=1:num_individuos
    if i <= num_individuos_select
        population(rank(i)).fitrel=fitness_1(i)/soma;
    else
        rank_excluir(i-num_individuos_select)=rank(i);

    end
end

evaluator = rank_excluir<best_gen;
best_gen = best_gen - sum(evaluator);

%eliminar ficheiros dos elementos ignorados
for k = 1:length(rank_excluir)
    nameDat = split(population(rank_excluir(k)).name,".");
    iname = nameDat(1) + "_input.txt";
    delete(population(rank_excluir(k)).name);
    delete(iname);
end

population(rank_excluir)=[];

end


function population = recomb(population,config,gen)
%Offspring
num_individuos = length(population);

fitrel=zeros(num_individuos,1);
a=(1:num_individuos); %a: indice do individuo
for i=1:num_individuos
    fitrel(i)= population(i).fitrel;
end
vect_aux=[a ; fitrel']; %tem os indices e a respectiva fiteness
%relativa (probabilidade de recombinar)

for i=1:config(3+config(21))
    n=num_individuos+i; %n da gera????o
    b{1}=randsrc(1,2,vect_aux); %extradorso: escolhe um 1x2 que tem os 2 indiv.
    b{2}=randsrc(1,2,vect_aux); %intradorso: escolhe um 1x2 que tem os 2 indiv.
    %     b{1}=randsrc(1,2,vect_aux);    Isto era para o caso em que eu
    %     recombino muito
    %     b{2}=randsrc(1,2,vect_aux);
    %     b{3}=randsrc(1,2,vect_aux);
    %     b{4}=randsrc(1,2,vect_aux);
    %Pontos fixos:
    %x
    offspring(1,1)=1;
    offspring(6,1)=0;
    offspring(11,1)=1;
    
    %y
    offspring(1,2)=0;
    offspring(6,2)=0;
    offspring(11,2)=0;
    
    %Recombina????o:
    for k=1:2
        %Extradorso
        offspring(2,k)=population(b{1}(1)).control(2,k)+0.002*randn;
        offspring(10,k)=population(b{1}(1)).control(10,k)+0.002*randn;
        
        offspring(3,k)=population(b{1}(2)).control(3,k)+0.002*randn;
        offspring(4,k)=population(b{1}(2)).control(4,k)+0.002*randn;
        
        %Intradorso
        offspring(5,k)=population(b{2}(1)).control(5,k)+0.002*randn;
        offspring(7,k)=population(b{2}(1)).control(7,k)+0.002*randn;
        
        offspring(8,k)=population(b{2}(2)).control(8,k)+0.002*randn;
        offspring(9,k)=population(b{2}(2)).control(9,k)+0.002*randn;
    end
    
    population(n).control = offspring;
    population(n).gen = gen;
    population(n).n = i;
    population(n).name = "createdAirfoils/perfil_" + int2str(gen) + "-" + int2str(i) + ".txt";
end
end
