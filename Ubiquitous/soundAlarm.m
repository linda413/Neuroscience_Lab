WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
Audio = audioplayer(WarnWave*20, 22050);
play(Audio);