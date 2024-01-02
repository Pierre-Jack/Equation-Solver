import math

xl = float(input("xl = "))
xu = float(input("xu = "))
lo = abs(xu - xl)
ea = 1E-2
fl = xl ** 2 - 9
fu = xu ** 2 - 9
significant_figures = 5
k = math.ceil(math.log2(lo / ea))
for i in range(0, k + 1):
    xr = (xl + xu) / float(2)
    xr = round(xr, -int(math.floor(math.log10(abs(xr)))) + (significant_figures - 1))
    fr = xr ** 2 - 9
    if fr < 0:
        xl = xr
    elif fr > 0:
        xu = xr
print("Bisection : x = " + str(xr))

for i in range(0, k + 1):
    xr = ((xl * fu) - (xu * fl)) / float(fu - fl)
    xr = round(xr, -int(math.floor(math.log10(abs(xr)))) + (significant_figures - 1))
    fr = xr ** 2 - 9
    if not fr == 0:
        round(fr, -int(math.floor(math.log10(abs(fr)))) + (significant_figures - 1))
    if fr < 0:
        xl = xr
        fl =  xl ** 2 - 9
        round(fl, -int(math.floor(math.log10(abs(fl)))) + (significant_figures - 1))
    elif fr > 0:
        xu = xr
        fu =  xu ** 2 - 9
        round(fu, -int(math.floor(math.log10(abs(fu)))) + (significant_figures - 1))
print("Regula-Falsi : x = " + str(xr))
