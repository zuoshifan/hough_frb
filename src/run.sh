#!/usr/bin/env bash

# generate low resolution sim and plot
python gen_sim_low_resolution.py -n -10 -m 10

# plot sim A
python plot_A_sim.py -n -200 -m 200
python plot_A3d_sim.py


# # plot image of pulsars
python plot_I_B0329.py -n 100 -m 220
python plot_I_B1929.py -n 100 --vmax 200
python plot_I_B2319.py -n 100 -m 225

# # plot image of FRBs
# python plot_I_010125.py
python plot_I_010125.py -n -0.4 -m 0.4
# python plot_I_010621.py
python plot_I_010621.py -n -0.4 -m 0.4
# python plot_I_010724.py
python plot_I_010724.py -n -1.1 -m 1.1
# python plot_I_110220.py
python plot_I_110220.py -n -0.4 -m 0.4
# python plot_I_110523.py
python plot_I_110523.py -n -9 -m 9
# python plot_I_110626.py
python plot_I_110626.py -n -0.4 -m 0.5
# python plot_I_110703.py
python plot_I_110703.py -n -0.4 -m 0.4
# python plot_I_120127.py
python plot_I_120127.py -n -0.4 -m 0.4
# python plot_I_140514.py
python plot_I_140514.py -n -10 -m 10

# plot accumulator of pulsars
python plot_A_B0329.py -n -240000 -m 240000
python plot_A_B1929.py -n -24500 -m 5000
python plot_A_B2319.py -n -80000 -m 80000

# plot accumulator of FRBs
# python plot_A_010125.py -n -1.6 -m 1.6
python plot_A_010125.py -n -3.6 -m 3.6
# python plot_A_010621.py -n -1.6 -m 1.6
python plot_A_010621.py -n -2.8 -m 2.8
# python plot_A_010724.py -n -40 -m 40
python plot_A_010724.py -n -60 -m 60
# python plot_A_110220.py -n -10 -m 10
python plot_A_110220.py -n -26 -m 26
python plot_A_110523.py -n -600 -m 600
python plot_A_110626.py -n -10 -m 10
python plot_A_110703.py -n -4 -m 4
plot_A_120127.py -n -4 -m 4
python plot_A_140514.py -n -160 -m 160

# plot PMR
python plot_pmr_sigmans_taus.py
python plot_pmr_taus_mus.py