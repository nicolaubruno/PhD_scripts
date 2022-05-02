# Libraries
import numpy as np
import matplotlib.pyplot as plt

# Parameters
F1 = 25.0 # mm
F2 = 25.0 # mm
d1 = 75.0 # mm
d2 = 80.4 # mm
d3 = 85.0 # mm
wavelength = 689e-6 # mm

# ABCD matrix 
# ---
def ABCD_matrix(F1, F2, d1, d2, d3):
	# Free propagation starting from mirror 2
	M1 = np.eye(2)
	M1[0][1] = d3

	# Reflection on the Mirror 1
	M2 = np.eye(2)
	M2[1][0] = -1 / F1

	# Free propagation starting from mirror 1
	M3 = np.eye(2)
	M3[0][1] = 2*d2 + d1

	# Reflection on the Mirror 2
	M4 = np.eye(2)
	M4[1][0] = -1 / F2

	# Result
	return np.dot(M4, np.dot(M3, np.dot(M2, M1)))

def get_q_prime_idx_from_radii(q_prime_arr, F1, F2):
	for i in range(q_prime_arr.shape[0]):
		if (q_prime_arr[i,0].real == F1) and (q_prime_arr[i,1].real == F2):
				break

	return i
# ---

# Available mirrors
F_arr = np.array([15])
F_arr = np.concatenate((F_arr, np.arange(25, 100, 25)))
F_arr = np.concatenate((F_arr, np.arange(100, 300, 50)))
F_arr = np.concatenate((F_arr, np.arange(500, 1250, 250)))

# Stability map
R_interval = np.arange(25, 1025, 25)
stability_map = np.zeros((R_interval.size, R_interval.size))
waist_arr = []

# Stability analisis considering the available mirrors
# ---
for i, F1 in enumerate(R_interval):
	for j, F2 in enumerate(R_interval):
		# Calculate ABCD Matrix
		M = ABCD_matrix(F1, F2, d1, d2, d3)

		# Check stability
		if abs(M[0][0] - M[1][1]) <= 2: 
			if (F1 in F_arr) and (F2 in F_arr): 
				# Set point in the stability map as waist
				q_prime = complex((M[1][1] - M[0][0]) / (2 * M[0][1]), - np.sqrt(4 - (M[0][0] - M[1][1])**2) / (2 * abs(M[0][1])))
				waist = -(wavelength * np.imag(q_prime)) / (np.pi * np.abs(q_prime)**2)*1e3 # micrometer
				stability_map[i][j] = waist
				waist_arr.append([F1, F2, waist])

# Waist array
waist_arr = np.array(waist_arr)
max_waist = np.max(waist_arr[:,2])
min_waist = np.min(waist_arr[:,2])
# ---

'''
for i, F1 in enumerate(R_interval):
	for j, F2 in enumerate(R_interval):
		if (F1 in np.real(q_prime_arr[:,0])) and (F2 in np.real(q_prime_arr[:,1])):
			q_i = get_q_prime_idx_from_radii(q_prime_arr, F1, F2)
			stability_map[i][j] = waist_arr[q_i] / max_waist
'''
# ---

# Plot stability graph
plt.clf()
plt.style.use("seaborn-white")
plt.rcParams.update({"font.size": 14, "figure.figsize": (18,5)})
fig, ax = plt.subplots(1,3)

# Subplot 01
map1 = ax[0].scatter(waist_arr[:, 0], waist_arr[:, 1], c=(waist_arr[:,2]), cmap="copper", marker="o", edgecolor="black", linewidth=1, s=50)

# Subplot 02
#--
map2 = ax[1].scatter(waist_arr[:, 0], waist_arr[:, 1], c=(waist_arr[:,2]), cmap="copper", marker="s", linewidth=0, s=450)

for waist in waist_arr:
	if (waist[0] > 0 and waist[0] < 275) and (waist[1] > 0 and waist[1] < 275):
		ax[1].text(waist[0], waist[1] - 1, int(waist[2]), color="white", ha="center", va="center")

ax[1].set_xlim(0, 275)
ax[1].set_ylim(0, 275)
# ---


# Subplot 03
# ---
min_F = 200
max_F = 1100
map3 = ax[2].scatter(waist_arr[:, 0], waist_arr[:, 1], c=(waist_arr[:,2]), cmap="copper", marker="s", linewidth=0, s=450)

for waist in waist_arr:
	if (waist[0] > min_F and waist[0] <max_F) and (waist[1] > min_F and waist[1] < max_F):
		ax[2].text(waist[0], waist[1] - 1, int(waist[2]), color="black", ha="center", va="center")

ax[2].set_xlim(min_F, max_F)
ax[2].set_ylim(min_F, max_F)
# ---

# Labels
cbar = plt.colorbar(map1)
cbar.set_label(r"$w_0\ [\mu m]$")

for i in range(3):
	ax[i].set_xlabel(r"$F_1\ [mm]$")
	ax[i].set_ylabel(r"$F_2\ [mm]$")
	ax[i].grid(linestyle="--")

plt.tight_layout()
plt.close(1)
plt.show()


