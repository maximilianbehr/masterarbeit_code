from src.instationary.chorin import Chorin


REs = [200]
refs = [3]
T = 10.0
savefreq = 10
dt = None

for ref in refs:
    for RE in REs:

        chorin = Chorin(ref, RE, T, savefreq, dt)
        chorin.solve()



