function value = get_NF_value(mycounter, max_frames)

%    sham_nf_value = randi(100,1, 1);
%    value = sham_nf_value;
    
    if mycounter <= round(max_frames/2) % before the middle point of the game
        value = 20;
    else
        value = 80;

    end
end

