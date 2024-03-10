%% Codigo que trabaja con la estación d10 para un nuevo análisis
clear all; close all; clc

% cargamos variables
load('datad10L2act.mat');
[latd10,lond10,altd10,name,time,velv,dirv,temp,hrel,humesp,pres,rglb] = load_est(datad10L2);
% Realizamos promedios horarios para las variables
N = length(time);
Nhoras = N/6;
Ndias = Nhoras/24;
velvh = nanmean(reshape(velv,6,Nhoras))';
temph = nanmean(reshape(temp,6,Nhoras))';
hrelh = nanmean(reshape(hrel,6,Nhoras))';
humesph = nanmean(reshape(humesp,6,Nhoras))';
presh = nanmean(reshape(pres,6,Nhoras))';
u = -velv.*sind(dirv);
v = -velv.*cosd(dirv);
uh = nanmean(reshape(u,6,Nhoras))';
vh = nanmean(reshape(v,6,Nhoras))';
dirvh = atan2d(uh,vh)+180;

% Ajustamos el tiempo
timehm = reshape(time,6,Nhoras)';
timeh = timehm(:,1);
timehd = reshape(timeh,24,Ndias)';
timed = timehd(:,1);

%% Trabajar los datos a niveles horarios, proyectando en la direccion de ref

% Direcciones a usar por inspeccion
dref1 = 43;    % noche
dref2 = 208;  % tarde

[Vref_noc,Vref_diur] = calc_vref(velv,dirv,dref1,dref2);

%% 1. Graficar series completas (horarias)

% desde donde parten los datos para el xlim
time1 = find(timeh == datenum(datetime(2013,1,1)));

% graficar velocidad
figure(1)
subplot(2,1,1)
plot(timeh,velvh,'b')
grid on
ylabel('V (m/s)')
title(['Serie de Tiempo de Velocidad de Viento a 20 m - Estación ',name])
ylim([0 25])
ticks = datenum(datetime(2013,1,1) + calmonths(0:12:72));
labels = datestr(ticks,'mmmyyyy');
set(gca, 'xtick', ticks, 'xtickLabel', labels);
xlim([timeh(time1) ticks(end)])
% graficar direccion
subplot(2,1,2)
plot(timeh,dirvh,'b.')
ylabel('Dir (°)')
title(['Serie de Tiempo de Dirección de Viento a 10 m - Estación ',name])
ylim([0 360])
yticks(0:90:360)
grid on
set(gca, 'xtick', ticks, 'xtickLabel', labels)
xlim([timeh(time1) ticks(end)])

% Graficar temperatura
figure(2)
subplot(2,1,1)
plot(timeh,temph,'b')
grid on
ylabel('T (°C)')
title(['Serie de Tiempo de Temperatura a 2 m - Estación ',name])
ylim([0 30])
set(gca, 'xtick', ticks, 'xtickLabel', labels);
xlim([timeh(time1) ticks(end)])
% graficar humedad esp
subplot(2,1,2)
plot(timeh,humesph,'b')
ylabel('q (g/kg)')
title(['Serie de Tiempo de Humedad Específica a 2 m - Estación ',name])
grid on
set(gca, 'xtick', ticks, 'xtickLabel', labels)
xlim([timeh(time1) ticks(end)])

cond2017_2meses = (month(timeh) >= 4 & month(timeh) <= 6) ...
    & (year(timeh) == 2017);
% cond2017_2meses = (month(timeh) >= 6 & month(timeh) <= 8) ...
%     & (year(timeh) == 2016);
timeh_2017_2meses = timeh(cond2017_2meses);
ticks2 = datenum(datetime(2017,4,1) + caldays(0:7:68));
% ticks2 = datenum(datetime(2016,6,1) + caldays(0:7:68));
labels2 = datestr(ticks2,'dd/mm');
velvh_2017_2meses = velvh(cond2017_2meses);
dirvh_2017_2meses = dirvh(cond2017_2meses);
temph_2017_2meses = temph(cond2017_2meses);
humesph_2017_2meses = humesph(cond2017_2meses);

% Serie de 2 meses
figure(3)
subplot(2,1,1)
hold on
plot(timeh_2017_2meses,velvh_2017_2meses,'b','linewidth',1)
grid on; ylabel('V (m/s)'); xlabel('Fechas')
% xlim([timeh_2017_4meses(1) timeh_2017_4meses(end)])
title('Serie de Tiempo de Velocidad para Abril-Mayo de 2017')
% title('Serie de Tiempo de Velocidad para Junio-Julio de 2016')
set(gca, 'xtick', ticks2, 'xtickLabel', labels2);
xlim([ticks2(1) ticks2(end)]); box on
subplot(2,1,2)
plot(timeh_2017_2meses,dirvh_2017_2meses,'.b','markersize',6)
grid on; ylabel('Dirección (°)'); xlabel('Fechas')
% xlim([timeh_2017_4meses(1) timeh_2017_4meses(end)])
title('Serie de Tiempo de Dirección de Viento para Abril-Mayo de 2017')
% title('Serie de Tiempo de Dirección de Viento para Junio-Julio de 2016')
set(gca, 'xtick', ticks2, 'xtickLabel', labels2);
xlim([ticks2(1) ticks2(end)])
yticks(0:90:360)
ylim([0 360])

figure(4)
subplot(2,1,1)
hold on
plot(timeh_2017_2meses,temph_2017_2meses,'b','linewidth',1)
grid on
ylabel('T (°C)')
title('Serie de Tiempo de Temperatura para Abril-Mayo de 2017')
ylim([0 30]); xlabel('Fechas')
set(gca, 'xtick', ticks2, 'xtickLabel', labels2);
% xlim([ticks2(1) ticks2(end)])
subplot(2,1,2)
plot(timeh_2017_2meses,humesph_2017_2meses,'b','linewidth',1)
grid on; ylabel('q (g/kg)')
xlim([timeh_2017_2meses(1) timeh_2017_2meses(end)])
title('Serie de Tiempo de Humedad Específica para Abril-Mayo de 2017')
set(gca, 'xtick', ticks2, 'xtickLabel', labels2);
% xlim([ticks2(1) ticks2(end)])

%% 2. Graficar proyeccion a nivel diario (global) nocturno y diurno

% desde donde empieza a contar
time2 = datenum(datetime(2013,1,1));%find(timed == datenum(datetime(2013,1,1)));
timed1 = timed(1:end-1);
Vref_noc = Vref_noc(2:end); %Vref_noc(2:end);
Vref_diur = Vref_diur(1:end-1);
Vref_dif = Vref_noc - Vref_diur;

%% Pequeño test para diagrama polar

uh_d10 = reshape(uh,24,Ndias)';
vh_d10 = reshape(vh,24,Ndias)';
unoc_d10 = nanmean(uh_d10(:,1:6),2);
vnoc_d10 = nanmean(vh_d10(:,1:6),2);
unoc_d10 = unoc_d10(1:end-1);
vnoc_d10 = vnoc_d10(1:end-1);
eventos = Vref_dif >= 10;
noeventos = Vref_dif <= 0;
mjja = (month(timed1) >= 5) & (month(timed1) <= 8); 
eventos_mjja = mjja & eventos;
noeventos_mjja = mjja & noeventos;
u_ev_d10 = unoc_d10(eventos_mjja);
v_ev_d10 = vnoc_d10(eventos_mjja);
u_noev_d10 = unoc_d10(noeventos_mjja);
v_noev_d10 = vnoc_d10(noeventos_mjja);

figure(20)
plot(u_ev_d10,v_ev_d10,'b.'); hold on
grid on; title('Diagrama Polar Estación d10 para Eventos en MJJA')
xlabel('u_{noc} (m/s)')
ylabel('v_{noc} (m/s)')

figure(21)
plot(u_noev_d10,v_noev_d10,'b.'); hold on
grid on; title('Diagrama Polar Estación d10 para Antieventos en MJJA')
xlabel('u_{noc} (m/s)')
ylabel('v_{noc} (m/s)')

%% Grafico de indices diurno y nocturno globales

figure(5)
subplot(2,1,1)
a1 = plot(timed1,Vref_noc,'b','linewidth',0.6); hold on
a2 = plot(timed1,Vref_diur,'r','linewidth',0.6);
plot(timed1,zeros(1,length(timed1)),'k')
grid on; legend([a1 a2],'Nocturno','Diurno'); ylabel('V_{ref} (m/s)')
title(['Serie de Tiempo de V_{ref} - Estación ',name])
xlim([time2 ticks(end)])
set(gca, 'xtick', ticks, 'xtickLabel', labels)
% Graficar la diferencia
subplot(2,1,2)
plot(timed1,Vref_dif,'b','linewidth',1); hold on
plot(timed1,zeros(1,length(timed1)),'k')
xlim([time2 ticks(end)])
grid on; ylabel('\Delta{V_{ref}} (m/s)')
title(['Serie de Tiempo de \Delta{V_{ref}} - Estación ',name])
set(gca, 'xtick', ticks, 'xtickLabel', labels)

% Temperatura por día, promedio diario
tempd = nanmean(reshape(temph,24,Ndias)',2); tempd = tempd(2:end);
humespd = nanmean(reshape(humesph,24,Ndias)',2); humespd = humespd(2:end);
presd = nanmean(reshape(presh,24,Ndias)',2); presd = presd(2:end);

% filtro = ~isnan(tempd);
inicio = find(timed == datenum(2013,1,24));
final = find(timed == datenum(2018,12,31));
filtro = inicio:final;
tempdn = tempd(filtro);
humespdn = humespd(filtro);
timed1n = timed1(filtro);
presdn = presd(filtro);
Vref_nocn = Vref_noc(filtro);
Vref_diurn = Vref_diur(filtro);

[tempdf] = Fourier(tempdn);
[humespdf] = Fourier(humespdn);
[presdf] = Fourier(presdn);
[Vref_nocf] = Fourier(Vref_nocn);
[Vref_diurf] = Fourier(Vref_diurn);
Vref_diff = Vref_nocf - Vref_diurf;

% Considerar que todas las series de aquí para abajo están con ciclo anual
% eliminado
figure(6)
subplot(3,1,1)
plot(timed1n,tempdf,'b','linewidth',1.1); hold on
grid on; ylabel('T(°C)')
title(['Serie de Tiempo de Temperatura Promedio Diario - Estación ',name])
set(gca, 'xtick', ticks, 'xtickLabel', labels)
xlim([time2 ticks(end)])
% Humedad específica promedio diario
subplot(3,1,2)
humespdf(humespdf < 0) = 0;
plot(timed1n,humespdf,'b','linewidth',1.1); hold on
grid on; ylabel('q (g/kg)')
title(['Serie de Tiempo de Humedad Específica Promedio Diario - Estación ',name])
set(gca, 'xtick', ticks, 'xtickLabel', labels)
xlim([time2 ticks(end)])
% Presión atmosférica promedio diario
subplot(3,1,3)
plot(timed1n,presdf,'b','linewidth',1.1); hold on
grid on; ylabel('p(hPa)')
title(['Serie de Tiempo de Presión Promedio Diario - Estación ',name])
set(gca, 'xtick', ticks, 'xtickLabel', labels)
xlim([time2 ticks(end)])

% Eliminar NaN para hacer el histograma
Vref_dif_snan = Vref_dif(~isnan(Vref_dif));
test = sum(Vref_dif_snan >= 10)/length(Vref_dif_snan);

figure(7)
edges = -14:2:36;
histogram(Vref_dif_snan,edges,'Normalization','Probability');
xlabel('\Delta{V_{ref}} (m/s)')
title(['Histograma Global de \Delta{V_{ref}} - Estación ', name])
ylabel('(%)'); xticks(-14:4:38)
yticks(0:0.02:0.2); ylim([0 0.2])
set(gca, 'YTick',yticks,'YTickLabel',yticks*100)
grid on; xlim([-14 34])

% Considerar el numero de veces que la velocidad supera el umbral de los
% 10 m/s en un mes. Obtener un porcentaje y hacer un boxplot
for i = 1:12
    vector = Vref_dif(month(timed1) == i);      
    largo = length(vector);
    if largo < 465                  % esto equipara el largo de los meses
        vector(largo:465) = NaN;
    end
    Vdif_box(:,i) = vector;        % nos servirá para un boxplot
    vector_sinan = vector(~isnan(vector));
    % porcentaje de dias que superan el umbral
    perc(i) = (sum(vector_sinan >= 10)/length(vector_sinan))*100;
end

for i = 2:2:12
    vect = Vref_dif(month(timed1) == i-1 | month(timed1) == i);
    largo = length(vect);
    if largo < 930                  % esto equipara el largo de los meses
        vect(largo:930) = NaN;
    end
    Vdif_box_bi(:,i/2) = vect;        % nos servirá para un boxplot
    vect_sinan = vect(~isnan(vect));
    perc_bi(i/2) = (sum(vect_sinan >= 10)/length(vect_sinan))*100;
end

for i = 3:3:12
    vec = Vref_dif((month(timed1) == i-2 | month(timed1) == i-1) ...
        | month(timed1) == i);
    largo = length(vec);
    if largo < 1380             % esto equipara el largo de los meses
        vec(largo:1380) = NaN;
    end
    Vdif_box_tri(:,i/3) = vec;        % nos servirá para un boxplot
    vec_sinan = vec(~isnan(vec));
    perc_tri(i/3) = (sum(vec_sinan <= 0)/length(vec_sinan))*100;
end

% Boxplot mensual
figure(8)
boxplot(Vdif_box)
ylabel('\Delta{V_{ref}} (m/s)')
xlabel('Meses')
title('Boxplot Global de \Delta{V_{ref}}')
grid on

% % Boxplot bimensual
% figure(9)
% boxplot(Vdif_box_bi)
% ylabel('\Delta{V_{ref}} (m/s)')
% meses = ["Ene-Feb","Mar-Abr","May-Jun","Jul-Ago","Sep-Oct","Nov-Dec"];
% xticklabels(meses)
% % horas = [24,06,12,18,24,06,12,18,24,06,12,18,24];
% % xticklabels(horas)
% title('Boxplot Global de \Delta{V_{ref}} Bimensual')
% grid on

% Boxplot trimensual
figure(9)
boxplot(Vdif_box_tri)
ylabel('\Delta{V_{ref}} (m/s)')
meses = ["DEF","MAM","JJA","SON"];
xticklabels(meses)
title('Boxplot Global de \Delta{V_{ref}} Trimestral')
grid on

% Hacer una regresión para la dispersión de ambos índices diurno y nocturno
X=[ones(length(Vref_diur),1) Vref_diur];
[b11,bint11,r11,rint11,stats11] = regress(Vref_noc,X); %r1 son los residuales

% Graficar dispersion entre ambos indices pero solo en MJJA
timevecd1 = datevec(timed1n);
cond_4meses = (timevecd1(:,1)>=2013 & timevecd1(:,1)<=2018) & ...
    (timevecd1(:,2)==5 | timevecd1(:,2)==6 | timevecd1(:,2)==7 ...
    | timevecd1(:,2)==8);
inicio1 = find(timed1 == datenum(2013,1,24));
fin1 = find(timed1 == datenum(2018,12,31));
Vref_noc2 = Vref_noc(inicio1:fin1);
Vref_diur2 = Vref_diur(inicio1:fin1);

Vref_noc4meses = Vref_noc2(cond_4meses);
Vref_diur4meses = Vref_diur2(cond_4meses);
temp_4meses = tempdf(cond_4meses);
humesp_4meses = humespdf(cond_4meses);
Vref_nocf4meses = Vref_nocf(cond_4meses);
Vref_diurf4meses = Vref_diurf(cond_4meses);
Vref_diff4meses = Vref_nocf4meses - Vref_diurf4meses;

X2 =[ones(length(Vref_diur4meses),1) Vref_diur4meses];
[b22,bint22,r22,rint22,stats22] = regress(Vref_noc4meses,X2);

% Graficar dispersion entre ambos indices, ver que tan anticorr están
figure(10)
% e1 = plot(Vref_diur,Vref_noc,'.b'); hold on
% e2 = plot(Vref_diur4meses,Vref_noc4meses,'.r');
size = 20; hold on;
e1 = scatter(Vref_diur,Vref_noc,size,'b');
e2 = scatter(Vref_diur4meses,Vref_noc4meses,size,'r');
plot(-20:15,zeros(1,length(-20:15)),'k')
plot(zeros(1,length(-5:20)),-5:20,'k')
plot(Vref_diur,Vref_diur*b11(2)+b11(1),'b','linewidth',1.1)
plot(Vref_diur4meses,Vref_diur4meses*b22(2)+b22(1),'r','linewidth',1.1)
axis equal
xlabel('V_{ref} diurno (m/s)'); ylabel('V_{ref} nocturno (m/s)'); grid on
title(['Dispersión de V_{ref} Diurno y Nocturno - Estación ',name])
legend(e2,'Meses MJJA')
correla = corrcoef(Vref_diur,Vref_noc,'Rows','pairwise');
% de forma normal, la correlacion es -0.63
% usando la noche siguiente, la correlación es de -0.67
% usando la noche previa, la correlación es de -0.42
[h,p] = ttest(Vref_diur,Vref_noc,'Alpha',0.05);

figure(11)
plot(Vref_diur4meses,Vref_noc4meses,'b.')
hold on
plot(-20:15,zeros(1,length(-20:15)),'k')
plot(zeros(1,length(-5:20)),-5:20,'k')
plot(Vref_diur4meses,Vref_diur4meses*b22(2)+b22(1),'r','linewidth',1.1)
xlabel('V_{ref} diurno (m/s)'); ylabel('V_{ref} nocturno (m/s)'); grid on
title(['Dispersión V_{ref} Diurno y Nocturno para MJJA - Estación ',name])
correla2 = corrcoef(Vref_diur4meses,Vref_noc4meses,'Rows','pairwise');

Vref_dif4meses = Vref_noc4meses - Vref_diur4meses;
events = Vref_dif4meses >= 10;
Vref_dif4meses_ev = Vref_dif4meses(events);

X3=[ones(length(temp_4meses),1) temp_4meses];
[b3,bint3,r3,rint3,stats3] = regress(Vref_dif4meses,X3);

figure(12)
scatter(temp_4meses,Vref_dif4meses,'b'); hold on
plot(temp_4meses,temp_4meses*b3(2)+b3(1),'r','linewidth',1.1)
plot(0:25,zeros(1,length(0:25)),'k')
ylabel('\Delta{V_{ref}} (m/s)'); xlabel('T (°C)'); grid on
title(['(a) Dispersión de \Delta{V_{ref}} y T en MJJA - Estación ',name])
correla3 = corrcoef(temp_4meses,Vref_dif4meses,'Rows','pairwise');
% corr entre indice noc y temp es -0.43
% corr entre indice diur y temp es 0.45
% corr entre indice dif y temp es -0.47 esta es la mayor7
T3 = table(temp_4meses,Vref_dif4meses);
mdl3 = fitlm(T3); %hace el modelo lineal
tbl3 = anova(mdl3); %realiza el anova

X4=[ones(length(humesp_4meses),1) humesp_4meses];
[b4,bint4,r4,rint4,stats4] = regress(Vref_dif4meses,X4);

figure(13)
scatter(humesp_4meses,Vref_dif4meses,'b'); hold on
plot(humesp_4meses,humesp_4meses*b4(2)+b4(1),'r','linewidth',1.1)
plot(0:6,zeros(1,length(0:6)),'k')
ylabel('\Delta{V_{ref}} (m/s)'); xlabel('q (g/kg)'); grid on
title(['(b) Dispersión de \Delta{V_{ref}} y q en MJJA - Estación ',name])
correla4 = corrcoef(humesp_4meses,Vref_dif4meses,'Rows','pairwise');
% corr entre indice noc y humesp es 0.30
% corr entre indice diur y humesp es -0.43
% corr entre indice dif y humesp es 0.41 esta es la mayor
T4 = table(humesp_4meses,Vref_dif4meses);
mdl4 = fitlm(T4); %hace el modelo lineal
tbl4 = anova(mdl4); %realiza el anova

X5=[ones(length(humesp_4meses),1) humesp_4meses];
[b5,bint5,r5,rint5,stats5] = regress(temp_4meses,X5);

figure(14)
scatter(humesp_4meses,temp_4meses,'b'); hold on
plot(humesp_4meses,humesp_4meses*b5(2)+b5(1),'r','linewidth',1.1)
xlabel('q (g/kg)'); ylabel('T (°C)'); grid on; ylim([0 22])
title(['(c) Dispersión de q y T en MJJA - Estación ',name])
correla5 = corrcoef(humesp_4meses,temp_4meses,'Rows','pairwise');
T5 = table(humesp_4meses,temp_4meses);
mdl5 = fitlm(T5); %hace el modelo lineal
tbl5 = anova(mdl5); %realiza el anova

%% 3. Graficar un año completo
% Año 2017

timevecd = datevec(timed);
cond_2017 = timevecd(:,1) == 2017;
Vref_noc2017 = Vref_noc(cond_2017);
Vref_diurno2017 = Vref_diur(cond_2017);
time_ini = datetime(2017,1,1);
time_fin = datetime(2017,12,31);
time2017 = time_ini:time_fin;

figure(15)
subplot(2,1,1)
b1 = plot(time2017,Vref_noc2017,'b','linewidth',1.1);
hold on
plot(time2017,Vref_noc2017,'.b','markersize',9)
b2 = plot(time2017,Vref_diurno2017,'r','linewidth',1.1);
plot(time2017,Vref_diurno2017,'.r','markersize',9)
plot(time2017,zeros(1,length(time2017)),'k')
set(gca,'fontsize',10)
legend([b1 b2],'Nocturno','Diurno')
ylabel('V_{ref} (m/s)','fontsize',10)
xlabel('Tiempo')
xlim([time_ini time_fin])
title(['Series de Tiempo de V_{ref} año 2017 - Estación ',name],'fontsize',11)
grid on
% Considerar la diferencia entre ambas series y graficarlas
subplot(2,1,2)
Vref_dif2017 = Vref_noc2017 - Vref_diurno2017;
plot(time2017,Vref_dif2017,'Color',[0.75 0 0.75],'linewidth',1.1); hold on
plot(time2017,Vref_dif2017,'.','Color',[0.75 0 0.75],'markersize',9)
plot(time2017,zeros(1,length(time2017)),'k')
set(gca,'fontsize',10); grid on
ylabel('\Delta{V_{ref}} (m/s)','fontsize',10)
xlabel('Tiempo')
title(['Serie de Tiempo de \Delta{V_{ref}} año 2017 - Estación ',name],'fontsize',11)
ylim([-10 35]); xlim([time_ini time_fin])

% Histograma de Delta Vref 2017
figure(16)
edges = -12:2:34;
histogram(Vref_dif2017(~isnan(Vref_dif2017)),edges,'Normalization','Probability')
xlabel('\Delta{V_{ref}} (m/s)')
title(['Histograma de \Delta{V_{ref}} año 2017 - ', name])
ylabel('(%)'); xticks(-12:4:34)
yticks(0:0.02:0.18); ylim([0 0.18])
set(gca, 'YTick',yticks,'YTickLabel',yticks*100)
grid on; xlim([-14 36])

% Ver que sucede en MJJA
cond_2017mjja = (timevecd1(:,1)==2017) & (timevecd1(:,2) >= 5 & ...
    timevecd1(:,2) <= 8);
Vref_noc2017mjja = Vref_noc2(cond_2017mjja);
Vref_diurno2017mjja = Vref_diur2(cond_2017mjja);
Vref_dif2017mjja = Vref_noc2017mjja - Vref_diurno2017mjja;
time_ini3 = datetime(2017,5,1);
time_fin3 = datetime(2017,8,31);
time2017mjja = time_ini3:time_fin3;

figure(17)
subplot(2,1,1)
b1 = plot(time2017mjja,Vref_noc2017mjja,'b','linewidth',1.1); hold on
plot(time2017mjja,Vref_noc2017mjja,'.b','markersize',10)
b2 = plot(time2017mjja,Vref_diurno2017mjja,'r','linewidth',1.1);
plot(time2017mjja,Vref_diurno2017mjja,'.r','markersize',10)
plot(time2017mjja,zeros(1,length(time2017mjja)),'k')
set(gca,'fontsize',10)
legend([b1 b2],'Nocturno','Diurno')
ylabel('V_{ref} (m/s)','fontsize',10)
xlabel('t (meses)')
xlim([time_ini3 time_fin3])
title(['Series de Tiempo para V_{ref} MJJA 2017 - Estación ',name], ...
    'fontsize',11)
grid on
% graficar la diferencia
subplot(2,1,2)
plot(time2017mjja,Vref_dif2017mjja,'b','linewidth',1.1); hold on
plot(time2017mjja,Vref_dif2017mjja,'.b','markersize',10)
set(gca,'fontsize',10)
ylabel('\Delta V_{ref} (m/s)','fontsize',10)
xlabel('t (meses)')
title(['Serie de Tiempo para \Delta V_{ref} MJJA 2017 - Estación ', ...
    name],'fontsize',11)
grid on; ylim([-10 35])
xlim([time_ini3 time_fin3])

% graficar temp y humedad
cond_2017amjja = (timevecd1(:,1)==2017) & (timevecd1(:,2) >= 4 & ...
    timevecd1(:,2) <= 8);
tempd_2017amjja = tempdf(cond_2017amjja);
humespd_2017amjja = humespdf(cond_2017amjja);
Vref_nocf2017amjja = Vref_nocf(cond_2017amjja);
Vref_diurf2017amjja = Vref_diurf(cond_2017amjja);
Vref_diff2017amjja = Vref_nocf2017amjja - Vref_diurf2017amjja;
time_ini4 = datetime(2017,4,1);
time_fin4 = datetime(2017,8,31);
time2017amjja = time_ini4:time_fin4;

figure(18)
subplot(2,1,1)
yyaxis left
c1 = plot(time2017amjja,tempd_2017amjja,'b','linewidth',1.1); hold on
plot(time2017amjja,tempd_2017amjja,'.b','markersize',9)
xlim([time_ini4 time_fin4])
set(gca,'fontsize',10); xlabel('Tiempo')
ylabel('T (°C)','fontsize',10); grid on
title(['Serie de Tiempo de Temperatura para Abril-Agosto de 2017 - Estación ',name], ...
    'fontsize',11)
yyaxis right
c2 = plot(time2017amjja,Vref_diff2017amjja,'r','linewidth',1.1); hold on
plot(time2017amjja,Vref_diff2017amjja,'.r','markersize',9)
ylabel('\Delta{V_{ref}} (m/s)')
% legend([c1 c2],'T','\Delta{V_{ref}}')
corr_temp = corrcoef(Vref_diff2017amjja,tempd_2017amjja,'Rows','pairwise');
subplot(2,1,2)
yyaxis left
d1 = plot(time2017amjja,humespd_2017amjja,'b','linewidth',1.1); hold on
plot(time2017amjja,humespd_2017amjja,'.b','markersize',9)
xlim([time_ini4 time_fin4])
set(gca,'fontsize',10); xlabel('Tiempo')
ylabel('q (g/kg)','fontsize',10); grid on
title(['Serie de Tiempo de Humedad Esp para Abril-Agosto de 2017 - Estación ',name], ...
    'fontsize',11)
yyaxis right
d2 = plot(time2017amjja,Vref_diff2017amjja,'r','linewidth',1.1); hold on
plot(time2017amjja,Vref_diff2017amjja,'.r','markersize',9)
ylabel('\Delta{V_{ref}} (m/s)')
% legend([d1 d2],'q','\Delta{V_{ref}}')
corr_hum = corrcoef(Vref_diff2017amjja,humespd_2017amjja,'Rows','pairwise');

%% 4. Graficar un año completo: 2015

% cond_2015 = timevecd(:,1)==2015;
% Vref_noc2015 = Vref_noc(cond_2015);
% Vref_diurno2015 = Vref_diur(cond_2015);
% time1s = datetime(2015,1,1);
% time2s = datetime(2015,12,31);
% time2015 = time1s:time2s;

% figure(19)
% subplot(2,1,1)
% c1 = plot(time2015,Vref_noc2015,'r','linewidth',1.1);
% hold on
% c2 = plot(time2015,Vref_diurno2015,'b','linewidth',1.1);
% plot(time2015,zeros(1,length(time2015)),'k')
% xlim([time1s time2s])
% legend([c1 c2],'Nocturno','Diurno')
% ylabel('V_{ref} (m/s)','fontsize',10)
% xlabel('t (meses)'); grid on
% title(['Series de Tiempo para V_{ref} año 2015 - Estación ',name],'fontsize',11)
% Vref_dif2015 = Vref_noc2015 - Vref_diurno2015;
% Vref_dif2015max = Vref_dif2015(Vref_dif2015 > 14);
% time2015max = time2015(Vref_dif2015 > 14);
% % Considerar la diferencia entre ambas series y graficarlas
% subplot(2,1,2)
% plot(time2015,Vref_dif2015,'b','linewidth',1.1); hold on
% plot(time2015,Vref_dif2015,'.b','markersize',10)
% %plot(time2015max,Vref_dif2015max,'r*')
% set(gca,'fontsize',10); grid on
% ylabel('\Delta V_{ref} (m/s)','fontsize',10)
% xlabel('t (meses)')
% title(['Serie de Tiempo para \Delta V_{ref} año 2015 - Estación ',name],'fontsize',11)
% xlim([time1s time2s]); ylim([-15 35])

% %% 5. Graficar promedio diario simple (nocturno)
% 
% Uref1 = -sind(dref1);
% Vref1 = -cosd(dref1);
% u = -velv.*sind(dirv);
% v = -velv.*cosd(dirv);
% velv_ref1 = u*Uref1 + v*Vref1;
% velv_refh1 = nanmean(reshape(velv_ref1,6,Nhoras))';
% velv_refhd1 = reshape(velv_refh1,24,Ndias)';
% Vref_noc_tot = nanmean(velv_refhd1,2);
% Vref_noc_tot = Vref_noc_tot(1:end-1);
% timestart = datenum(datetime(2013,1,1));
% 
% Vref_noctotmax = Vref_noc_tot(Vref_noc_tot>8);
% timedmax = timed(Vref_noc_tot>=8);
% figure()
% % promedio indice nocturno
% subplot(2,1,1)
% plot(timed1,Vref_noc_tot,'b','linewidth',1.1);
% hold on
% plot(timed1,zeros(1,length(timed1)),'k')
% %plot(timedmax,Vref_noctotmax,'r*')
% grid on
% ylabel('V_{ref} tot (m/s)')
% title(['Serie Completa de V_{ref} Total - Estación ',name])
% set(gca, 'xtick', ticks, 'xtickLabel', labels)
% xlim([timed1(time2) ticks(end)])
% 
% figure()
% Vref_noc_tot_snan = Vref_noc_tot(~isnan(Vref_noc_tot));
% histogram(Vref_noc_tot_snan,20,'Normalization','probability')
% xlabel('V_{ref} tot noc (m/s)')
% title(['Histograma de V_{ref} Total Global - ', name])
% grid on
% ylabel('(%)')
% yticks(0:0.02:0.16)
% ylim([0 0.16])
% set(gca, 'YTick',yticks,'YTickLabel',yticks*100)
% 
% for i = 1:12
%     vec = Vref_noc_tot(month(timed1) == i);
%     largo = length(vec);
%     if largo < 465                  % esto equipara el largo de los meses
%         vec(largo:465) = NaN;
%     end
%     Vref_noctot_box(:,i) = vec;        % nos servirá para un boxplot
%     vec_sinan2 = vec(~isnan(vec));
%     perc2(i) = (sum(vec_sinan2 >= 8)/length(vec_sinan2))*100;
%     % porcentaje de dias que superan el umbral
% end
% 
% % Boxplot
% figure()
% boxplot(Vref_noctot_box)
% ylabel('V_{ref} tot (m/s)')
% xlabel('Meses')
% title('Boxplot Global de V_{ref} Total')
% grid on
% ylim([-7 20])

% %% 6. Graficar series de 5 meses: mayo, junio, julio, agosto y septiembre de 2016 para Vref
% 
% timevech = datevec(timeh);
% cond_2016mes5 = (timevecd(:,1)==2016) & (timevecd(:,2)==5 | timevecd(:,2)==6 ...
%     | timevecd(:,2)==7 | timevecd(:,2)==8 | timevecd(:,2)==9);
% vref_2016noc = Vref_noc(cond_2016mes5);
% vref_2016diur = Vref_diur(cond_2016mes5);
% time_ini2 = datetime(2016,5,1);
% time_fin2 = datetime(2016,9,30);
% time_2016 = time_ini2:time_fin2;
% 
% figure()
% % indices nocturno y diurno
% subplot(2,1,1)
% d1 = plot(time_2016,vref_2016noc,'b','linewidth',1.1);
% hold on
% d2 = plot(time_2016,vref_2016diur,'r','linewidth',1.1);
% plot(time_2016,zeros(1,length(vref_2016noc)),'k')
% xlim([time_2016(1) time_2016(end)])
% ylabel('V_{ref} (m/s)')
% title(['Serie de V_{ref} (MJJAS de 2016) - Estación ',name])
% %dateFormat = 19;
% %datetick('x',dateFormat,'keepticks')
% %set(ax,'fontsize',9)
% grid on
% legend([d1 d2],'Nocturno','Diurno')
% % con diferencias
% vref_2016dif = vref_2016noc - vref_2016diur;
% subplot(2,1,2)
% plot(time_2016,vref_2016dif,'b','linewidth',1.1)
% hold on
% plot(time_2016,zeros(1,length(vref_2016noc)),'k')
% xlim([time_2016(1) time_2016(end)])
% grid on
% ylabel('\Delta{V_{ref}} (m/s)')
% title(['Serie de \Delta{V_{ref}} (MJJAS de 2016) - Estación ',name])
% %set(ax,'fontsize',9)
% %datetick('x',dateFormat,'keepticks')

% figure()
% edges = -7:1:35;
% histogram(vref_2016dif,edges,'Normalization','probability')
% xlabel('\Delta{V_{ref}} (m/s)')
% title(['Histograma de \Delta{V_{ref}} (MJJAS de 2016) - ', name])
% ylabel('(%)')
% xticks(-10:3:38)
% xlim([-10 38])
% yticks(0:0.05:0.20)
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
% ylim([0 0.2])
% grid on
% 
% % abril = vref_2017dif(1:30);
% % mayo = vref_2017dif(31:61);
% % junio = vref_2017dif(62:91);
% % perc_abril = (sum(abril > 14)/30)*100;
% % perc_mayo = (sum(mayo > 14)/31)*100;
% % perc_junio = (sum(junio > 14)/30)*100;

