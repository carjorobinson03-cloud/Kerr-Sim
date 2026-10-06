import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

d = np.load("render_linear.npz")
lum, r_cam, fov = d["lum"], float(d["r_cam"]), float(d["fov"])
H = lum.shape[0]

M_per_px = r_cam / np.sqrt(1 - 2/r_cam) * 2*np.tan(np.radians(fov/2)) / H #convert FOV deg into uas (microarcseconds)
uas_per_px = M_per_px * 3.8 # theta_g for M87* EHT Paper VI
sigma_px = (20 / 2.355) / uas_per_px  # 20 uas FWHM to Gaussian sigma in pixels
print(f"{M_per_px:.4f} M/px, {uas_per_px:.3f} uas/px, "
      f"20 uas FWHM = {20/uas_per_px:.0f} px (sigma {sigma_px:.0f} px)")

I = lum[::-1]                          
s = min(I.shape); y0 = (I.shape[0]-s)//2; x0 = (I.shape[1]-s)//2
I = I[y0:y0+s, x0:x0+s]
I = I / I.sum()
Ib = gaussian_filter(I, sigma_px, mode="constant")

half = s * uas_per_px / 2
ext = [-half, half, -half, half]
fig, ax = plt.subplots(1, 2, figsize=(11, 5))
ax[0].imshow(np.sqrt(I), extent=ext, cmap="afmhot"); ax[0].set_title("render (linear, sqrt scale)")
ax[1].imshow(Ib, extent=ext, cmap="afmhot");         ax[1].set_title("blurred to 20 uas FWHM")
for a in ax: a.set_xlabel("uas"); a.set_ylabel("uas")
plt.tight_layout(); plt.show()