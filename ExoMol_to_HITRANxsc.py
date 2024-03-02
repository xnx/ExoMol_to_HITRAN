"""
A script to convert ExoMol-style cross sections (those
produced by the ExoMol website's cross section service)
to HITRAN .xsc format.

"""

import os
import sys

sigma_path = sys.argv[1]
sigma_dir = os.path.dirname(sigma_path)
sigma_basename = os.path.basename(sigma_path)
sigma_filestem = os.path.splitext(sigma_basename)[0]
iso_formula, nurange, T, s_dnu = sigma_filestem.split('_')
s_numin, s_numax = nurange.split('-')
p = 0

numin, numax, dnu = map(float, [s_numin, s_numax, s_dnu])
if T[-1] == 'K':
    T = T[:-1]

npts_expected = round((numax - numin)/dnu) + 1

print(f"Isotopologue: {iso_formula}")
print(f"numin = {numin} cm-1")
print(f"numax = {numax} cm-1")
print(f"dnu = {dnu} cm-1")
print(f"Expected number of data points = {npts_expected}")
print(f"Temperature = {T} K")
print(f"Pressure = {p} Torr (no collisional broadening)")

# Read in the .sigma data, first determining if the file
# lines are nu, sigma or just sigma on its own.
with open(sigma_path) as fi:
    line = fi.readline().strip()
    fields = line.split()
    # Index at 1 for two-column layout or 0 for one-column
    idx = len(fields) == 2
    sigma = [float(fields[idx])]
    for line in fi.readlines():
        sigma.append(float(line.split()[idx]))

assert len(sigma) == (npts := npts_expected)
sigma_max = max(sigma)

xsc_name = f"{iso_formula}_{T}_{p}_{s_numin}-{s_numax}_XX.xsc"

xsc_path = os.path.join(sigma_dir, xsc_name)
print("Writing HITRAN-formatted cross section file to:")
print(xsc_path)

nrows = npts // 10
nextra = npts % 10
# If the number of points isn't a multiple of 10 then
# pad the final row with the right number of zeros until
# it is.
if nextra:
    nrows += 1
    sigma += [0.] * (10 - nextra)

with open(xsc_path, 'w') as fo:
    print(f"{iso_formula:>20s}{s_numin:>10s}{s_numax:>10s}{npts:>7d}{T:>7s}{p:>6.1f}{sigma_max:>10.3E}", file=fo)
    for i in range(nrows):
        print(''.join([f'{y:>10.3E}' for y in sigma[10*i:(10*i+10)]]), file=fo)
