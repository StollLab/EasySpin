function ok = test()

h = planck;
href = 6.62607015e-34; % CODATA 2018 / SI 2019 value; exact

ok = areequal(h,href,1e-9,'rel');
