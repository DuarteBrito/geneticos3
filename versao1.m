% perfil inicial
perfil=importdata('naca0009.txt', ' ');
x=perfil(:,1);
y=perfil(:,2);
extra = perfil(1:fix(length(perfil)/2),:);
intra = perfil(fix(length(perfil)/2)+1:length(perfil),:);
 
figure()
hold on
plot (x,y)


% a's iniciais
t1s = linspace(0.1,0.9999999,11);
t2 = 5;
as1 = rand(11,1)/1000;
as2 = -rand(11,1)/1000;
x_ = flip(extra(:,1));
%figure
for cnt = 1:length(t1s)
    t1 = t1s(cnt);
    a1 = as1(cnt);
    %x_ = flip(extra(:,1));
    %y_ = x_;
    y_ = a1 * sin(pi*x_.^(log(0.5)/log(t1))).^t2;
    plot(x_,y_)
    hold on
    extra(:,2) = extra(:,2) + y_;
end
x_ = flip(intra(:,1));
for cnt = 1:length(t1s)
    t1 = t1s(cnt);
    a2 = as2(cnt);
    %x_ = flip(extra(:,1));
    %y_ = x_;
    y_ = a2 * sin(pi*x_.^(log(0.5)/log(t1))).^t2;
    plot(x_,y_)
    hold on
    intra(:,2) = intra(:,2) + y_;
end
perfil_ = [extra ; intra];
x=perfil_(:,1);
y=perfil_(:,2);
plot (x,y)

