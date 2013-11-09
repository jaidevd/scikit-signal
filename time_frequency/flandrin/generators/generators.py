import numpy as np



def fmlin(N, fnormi=0.0, fnormf=0.5, t0=None):

    if t0 is None:
        t0 = np.round(N/2)

    if N <= 0:
        raise TypeError

    elif (np.abs(fnormi) > 0.5) or (np.abs(fnormf) > 0.5):
        raise TypeError

    else:
        y = np.arange(N);
        y = fnormi * (y - t0) + ((fnormf - fnormi)/(2.0*(N - 1))) * \
            ((y - 1)**2 - (t0 - 1)**2)
        y = np.exp(1j*2.0*np.pi*y)
        y = y/(t0)
        iflaw = np.linspace(fnormi, fnormf, N)

        return y, iflaw


def amgauss(N, t0=None, T=None):
    if t0 is None:
        t0 = np.round(N/2)

    if T is None:
        T = 2*np.sqrt(N);


    if N <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(N) - t0;
        y = np.exp(-((tmt0/T)**2)*np.pi)
        return y


if __name__ == "__main__":
    from matplotlib.pyplot import plot, show
    y = amgauss(256)
    plot(y)
    show()
