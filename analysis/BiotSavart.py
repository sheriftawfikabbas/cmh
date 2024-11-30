import numpy as np
C = 1.6e-19
bond = 3e-10
I = 1
mu0 = 4*np.pi*1e-7

class BiotSavart:

    def H(self, p0, z0, h):
        Hx = p0*np.log(1+(h**2+2*z0*h)/(p0**2+z0**2))+2*(z0+h) * \
            np.arctan(p0/(z0+h))-2*z0*np.arctan(p0/z0)
        # Hx = Hx/2
        Hz = z0*np.log(1+p0**2/z0**2)-(z0+h)*np.log(1+p0**2/(z0+h)
                                                    ** 2)-2*p0*(np.arctan((z0+h)/p0)-np.arctan(z0/p0))
        # Hz = Hz/2
        return Hx, Hz

    def B(self, x, z, w1, w2, h):
        A = (abs(w1)+abs(w2))*(abs(h))
        BS_X = I/A*mu0/(4*np.pi)*(self.H(w1+x, z, h)[0]-self.H(x+w2, z, h)[0])
        BS_Z = I/A*mu0/(4*np.pi)*(self.H(w1+x, z, h)[1]-self.H(x+w2, z, h)[1])
        return BS_X, BS_Z
