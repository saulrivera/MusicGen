%%% Algoritmo Genetico Musical
clc; clear; close all;

%% Variables de entrada

not = 30;       %numero de notas
ton = 20;       %numero de tonos

pc = 0.75;      %probabilidad de cruce
pm = 0.1;       %probabilidad de mutacion

maxg = 100;     % maxgimas generaciones
tmn = 100;     %tiempo maximo del tono
tmt = 200;        %tiempo maximo notas

%% Algoritmo Genetico

% mt = randi(12, not, ton);       %notas aleatorias
% t = round(rand(tmt, ton));     %tiempos aleatorios  
% sil = randi(tmn, not, ton);       %duracion de notas aleatorias

load('matrizp.mat')     % Cargar workspace pasados

notas_sonidos = @(frq,dur) cos(2*pi* (1:dur)/2048 * (447*2.^((frq-1)/12)));
nota = {'A' 'A#' 'B' 'C' 'C#' 'D' 'D#' 'E' 'F' 'F#' 'G' 'G#'}; 

fitp = zeros(1, 10);    %declarando vector de costo de notas
fitp2 = zeros(1, 10);   %declarando vector de costo de velocidad
fitp3 = zeros(1, 10);   %declarando vector de costo de duracion de notas

for p = 1:maxg
    
    clc
    fprintf(['Generación número ', num2str(p), '.'])    %imprime en la pantalla el numero de generacion

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Obtencion del peso
    
    cn = nota(mt);  %transforma el vector mt en notas

    for o = 1:ton
        cancion = cn(:, o)';
        for ind = 1:length(cancion)
            indx = strcmp(cancion(ind), nota);
            songidx(ind) = find(indx); %#ok<SAGROW>
        end    

        dur = sum(t(:, o))*2048/100;
        toca_cancion = [];

        for ind = 1:length(songidx)
           toca_cancion = [toca_cancion; [notas_sonidos(songidx(ind), dur)  zeros(1,sil(ind, o))]']; %#ok<AGROW>
        end
        
        fprintf(['\nTono ', num2str(o)])    %imprime el numero de tono que se esta reproduciendo
        tic
        soundsc(toca_cancion, 2048)
        
        v = 1;
        
        while v == 1
            v2 = input('\n¿Deseas repetir el tono? (Si = 1, No = 0): ');
            if v2 == 1
                soundsc(toca_cancion, 2048)
            elseif v2 == 0
                v = 0;
            end
        end
toc
        fitp(o) = input('¿Qué tal las notas? (0 - 10): ');
        fitp3(o) = input('¿Qué tal el ritmo? (0 - 10): ');
        fitp2(o) = input('¿Qué tal la velocidad? (0 - 10): ');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Seleccion Natural

    NP0 = mt(:, fitp > mean(fitp));     %obtiene solo aquellas notas mayores al promedio
    NQ0 = t(:, fitp2 > mean(fitp2));     %obtiene solo aquellas velocidades mayores al promedio
    NR0 = sil(:, fitp3 > mean(fitp3));   %obtiene solo aquellas duraciones mayores al promedio
    
    fitp = fitp(fitp > mean(fitp));             %actualizan los pesos
    fitp2 = fitp2(fitp2 > mean(fitp2));
    fitp3 = fitp3(fitp3 > mean(fitp3));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cruze de especies
    
    NP1 = zeros(size(mt));
    NQ1 = zeros(size(t));
    NR1 = zeros(size(sil));
    
    ind2 = randi(size(NP0, 2), 2, ton);
    ind3 = randi(size(NQ0, 2), 2, ton);
    ind4 = randi(size(NR0, 2), 2, ton);
    
    for o = 1:2:ton
        
        n1 = ind2(find(fitp(ind2(:, o)) == max(fitp(ind2(:, o))), 1), o);
        n2 = ind2(find(fitp(ind2(:, o + 1)) == max(fitp(ind2(:, o + 1))), 1), o + 1);
        
        n3 = ind3(find(fitp2(ind3(:, o)) == max(fitp2(ind3(:, o))), 1), o);
        n4 = ind3(find(fitp2(ind3(:, o + 1)) == max(fitp2(ind3(:, o + 1))), 1), o + 1);
        
        n5 = ind4(find(fitp3(ind4(:, o)) == max(fitp3(ind4(:, o))), 1), o);
        n6 = ind4(find(fitp3(ind4(:, o + 1)) == max(fitp3(ind4(:, o + 1))), 1), o + 1);

        g = randi(not);         % punto de corte notas
        g2 = randi(tmt);        %punto de corte tiempo tono
        g3 = randi(not);        % punto de corte tiempo entre notas
        
        NP1(:, o) = [NP0(1:g, n1); NP0(g + 1:end, n2)];
        NP1(:, o + 1) = [NP0(1:g, n2); NP0(g + 1:end, n1)];
        
        NQ1(:, o) = [NQ0(1:g2, n3); NQ0(g2 + 1:end, n4)];
        NQ1(:, o + 1) = [NQ0(1:g2, n4); NQ0(g2 + 1:end, n3)];
        
        NR1(:, o) = [NR0(1:g3, n5); NR0(g3 + 1:end, n6)];
        NR1(:, o + 1) = [NR0(1:g3, n6); NR0(g3 + 1:end, n5)];
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mutacion
    
    for o = 1:ton
        if rand <= pm
            l = randi(not);
            for u = 1:l
                NP1(randi(not), o) = randi(12);
                NQ1(randi(tmt), o) = randi(2) - 1;
                NR1(randi(not), o) = randi(tmn);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Elitismo
    
    [~, ind] = sort(fitp, 'descend');
    [~, ind22] = sort(fitp, 'descend');
    [~, ind33] = sort(fitp, 'descend');
    NP1(:, 1) = mt(:, ind(1));
    NQ1(:, 1) = t(:, ind22(1));
    NR1(:, 1) = sil(:, ind33(1));

    mt = NP1;
    t = NQ1;
    sil = NR1;
    
    save('matrizp.mat', 'mt', 't', 'sil');
    
end