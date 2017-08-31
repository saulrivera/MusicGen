function crearaudio(matriz, idsonido)

    load(matriz);

    notas_sonidos = @(frq, dur) cos(2*pi* (1:dur)/2048 * (447*2.^((frq-1)/12)));

    dur = sum(t(:, idsonido))*2048/100;

    songidx = mt(:, idsonido);
    
    toca_cancion = [];

    for ind = 1:length(songidx)
    toca_cancion = [toca_cancion; [notas_sonidos(songidx(ind), dur)  zeros(1,sil(ind, idsonido))]']; %#ok<AGROW>
    end
    
    name = ['audio', num2str(idsonido), '.wav'];
    
    audiowrite(name, toca_cancion, 2048)
    
    so = ['Audio creado con el nombre: ', name, '\n\n'];
    
    fprintf(so)
    
end