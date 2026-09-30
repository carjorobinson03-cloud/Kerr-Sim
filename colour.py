import numpy as np
import os
from scipy.integrate import cumulative_trapezoid


if not os.path.exists("bbtorgb.bin"):
    print("Baking Color Binary (One time)...")

def Flux_funcr(r_n, M_val, a_val, scale):
    """
    Takes nodular values for flux at points across disc, radius of nodular position is given by r_n.
    This will power the boltzmann relation to temperature and thus colour to be utilized in a lookup table.
    """
    #emission profile bit. ***THESE ARE ALL PROGRADE TAKING UPPER SIGNS; BPT73
    #Omega
    omeg_kepl = (np.sqrt(M_val)) / (r_n**(3/2) + a_val*np.sqrt(M_val))
    #partial Omega
    partial_r_Omeg = (-3.0*np.sqrt(M_val)*np.sqrt(r_n)) / (2.0*(r_n**(3.0/2.0)+a_val*np.sqrt(M_val))**2)
    #L dagger
    num_ldag = (np.sqrt(M_val)*(r_n**2 - 2.0*a_val*np.sqrt(M_val)*np.sqrt(r_n) + a_val**2))
    denom_ldag = (r_n**(3/4)*np.sqrt(r_n**(3/2) - 3.0*M_val*np.sqrt(r_n) + 2.0*a_val*np.sqrt(M_val)))
    Ldagger = num_ldag / denom_ldag
    #partial r Ldagger
    n_Ldagr_1 = np.sqrt(M_val)*r_n**4 - 6.0*M_val**(3/2)*r_n**3 + 8.0*M_val*a_val*r_n**(5/2) - 3.0*np.sqrt(M_val)*a_val**2*r_n**2
    n_Ldagr_2 = r_n*(M_val*a_val*r_n**(3/2) - 3.0*M_val**2*a_val*np.sqrt(r_n) + 8.0*M_val**(3/2)*a_val**2)
    n_Ldagr_3 = -3.0*M_val**2*a_val*r_n**(3/2) - 3.0*M_val*a_val**3*np.sqrt(r_n)
    numerator_Ldagr = n_Ldagr_1 + n_Ldagr_2 + n_Ldagr_3
    denominator_Ldagr = 2.0*r_n**(9/4)*((r_n**(3/2)-3.0*M_val*np.sqrt(r_n)+2.0*np.sqrt(M_val)*a_val)**(3/2))
    Ldag_partial_r = numerator_Ldagr / denominator_Ldagr
    # E dagger
    num_edag = (r_n**(3/2) - 2.0*M_val*np.sqrt(r_n) + a_val*np.sqrt(M_val))
    denom_edag = (r_n**(3/4)*np.sqrt(r_n**(3/2) - 3.0*M_val*np.sqrt(r_n) + 2.0*a_val*np.sqrt(M_val)))
    Edag = num_edag / denom_edag

    H = (Edag - omeg_kepl*Ldagger) * Ldag_partial_r #integrand
    I = cumulative_trapezoid(H, r_n, initial=0.0) # trapezoidal integration method because this would be a disaster to do regularily
    pref = -partial_r_Omeg / ((Edag - omeg_kepl*Ldagger)**2)
    return scale * (1.0/r_n) * pref * I     # sqrt(-g) = r, simple sympy stuff

#WSS13 and #CIE1931, these are experimentally accepted rgb coding approaches

def g(x, mu, s1, s2):
    S = np.where(x < mu, s1, s2)
    return np.exp(-0.5*((x - mu)/ S)**2) 

def cmf(lam):                              
    x = 1.056*g(lam,599.8,37.9,31.0) + 0.362*g(lam,442.0,16.0,26.7) - 0.065*g(lam,501.1,20.4,26.2)
    y = 0.821*g(lam,568.8,46.9,40.5) + 0.286*g(lam,530.9,16.3,31.1)
    z = 1.217*g(lam,437.0,11.8,36.0) + 0.681*g(lam,459.0,26.0,13.8)
    return x, y, z

lamda = np.linspace(380, 780, 401) #visible light spectrum, 401 nodes
lamda_m = lamda * 1e-9 #convert to meters
h, c, kB = 6.626e-34, 2.998e8, 1.381e-23
xb, yb, zb = cmf(lamda)

sRGB = np.array([[ 3.2406,-1.5372,-0.4986],
                [-0.9689, 1.8758, 0.0415],
                [ 0.0557,-0.2040, 1.0570]]) #standard d65 sRGB


def bb_to_rgb(T):
    #Plank's Law
    B = 1.0 / (lamda_m**5 * (np.exp(h*c/(lamda_m*kB*T)) - 1.0))
    #RGB code from integral, trapezoidal again
    X = np.trapezoid(B*xb, lamda)
    Y = np.trapezoid(B*yb, lamda) 
    Z = np.trapezoid(B*zb, lamda)
    xyz = np.array([X, Y, Z]) / Y
    #XYZ -> sRGB
    rgb = sRGB @ xyz
    rgb = np.clip(rgb, 0, None) #kill negatives
    rgb = rgb / max(rgb.max(), 1e-9) #normalize for brightest channel
    return rgb


T_LUT_MIN = 1000.0
T_LUT_MAX = 20000.0
T_LUT_N = 1024
T_grid = np.linspace(T_LUT_MIN, T_LUT_MAX, T_LUT_N)
blackbodtoRGB = np.array([bb_to_rgb(T) for T in T_grid], dtype=np.float32)

blackbodtoRGB.tofile("bbtorgb.bin")
print("colour baked")