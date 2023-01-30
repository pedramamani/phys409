from dynamics_config import *


def dynamics(I):
    Nk = np.zeros(shape=(N, K_LEVELS))
    Ni = np.zeros(shape=(N, I_LEVELS))
    Abs = np.zeros(shape=N)

    Pik, Gik = PIK_OVER_I * I, GIK
    Gi = np.sum(Gik, axis=1)
    Pi = np.sum(Pik, axis=1)
    Pk = np.sum(Pik, axis=0)

    Nk[0] = [1 / K_LEVELS] * K_LEVELS
    Abs0 = np.dot(Pk, Nk[0])
    Abs[0] = 1

    for n in range(N-1):
        Nk[n + 1] = Nk[n] + (Ni[n] @ Pik - np.multiply(Pk, Nk[n]) + Ni[n] @ Gik) * DT
        Ni[n + 1] = Ni[n] + (Nk[n] @ Pik.T - np.multiply(Pi, Ni[n]) - np.multiply(Gi, Ni[n])) * DT
        Abs[n + 1] = np.dot(Pk, Nk[n + 1]) / Abs0

    return np.linspace(0, T, N), Nk.T, Ni.T, Abs

