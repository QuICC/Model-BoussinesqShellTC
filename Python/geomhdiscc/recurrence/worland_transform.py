import sympy as sy

n = sy.Symbol('N')
a = sy.Symbol('a')
b = sy.Symbol('b')

#print("Natural normalization of Worland polynomials")
#print("===========================================")
# General Pn
Pnc0 = -(n+a-1)*(n+b-1)*(2*n+a+b)/(n*(n+a+b)*(2*n+a+b-2))
Pnc1 = (2*n+a+b-1)*(2*n+a+b)/(2*n*(n+a+b))
Pnc2 = (a*a - b*b)*(2*n+a+b-1)/(2*n*(n+a+b)*(2*n+a+b-2))
Pnc3 = 1
# P1
P1c0 = (a+b+2)/2
P1c1 = (a-b)/2
P1c2 = sy.Rational(1,2)
# P0
P0c0 = 1
#print("Polynomial recurrence:")
#print((Pnc0/Pnc3).simplify().factor())
#print((Pnc1/Pnc3).simplify().factor())
#print((Pnc2/Pnc3).simplify().factor())
#print(Pnc3)
#print("Polynomial n = 1 recurrence:")
#print((P1c0/P1c2).simplify().factor())
#print((P1c1/P1c2).simplify().factor())
#print(P1c2.simplify().factor())
#print("Polynomial n = 0 recurrence:")
#print(P0c0)
#print("-------------------------------------------")
# General DPn
DPnc0 = (n+a+b)*Pnc0/(n+a+b-2)
DPnc1 = (n+a+b)*Pnc1/(n+a+b-1)
DPnc2 = (n+a+b)*Pnc2/(n+a+b-1)
DPnc3 = 1/(n+a+b-1)
# DP1
DP1c0 = (a+b+1)*P1c0/(a+b)
DP1c1 = (a+b+1)*P1c1/(a+b)
DP1c2 = (a+b+1)/(2*(a+b))
# DP0
DP0c0 = 2*(a+b)
#print("First derivative recurrence:")
#print((DPnc0/DPnc3).simplify().factor())
#print((DPnc1/DPnc3).simplify().factor())
#print((DPnc2/DPnc3).simplify().factor())
#print(DPnc3)
#print("First derivative n = 1 recurrence:")
#print((DP1c0/DP1c2).simplify().factor())
#print((DP1c1/DP1c2).simplify().factor())
#print(DP1c2.simplify().factor())
#print("First derivative n = 0 recurrence:")
#print(DP0c0.simplify().factor())
#print("-------------------------------------------")
## General D2Pn
#D2Pnc0 = (n+a+b)*(n+a+b-1)*Pnc0/((n+a+b-2)*(n+a+b-3))
#D2Pnc1 = (n+a+b)*(n+a+b-1)*Pnc1/((n+a+b-1)*(n+a+b-2))
#D2Pnc2 = (n+a+b)*(n+a+b-1)*Pnc2/((n+a+b-1)*(n+a+b-2))
#D2Pnc3 = 1/(n+a+b-2)
## P1
#D2P1c0 = (a+b+1)*(a+b)*P1c0/((a+b)*(a+b-1))
#D2P1c1 = (a+b+1)*(a+b)*P1c1/((a+b)*(a+b-1))
#D2P1c2 = (a+b+1)/(2*(a+b-1))
## P0
#D2P0c0 = 4*(a+b)*(a+b-1)
#print("Second derivative recurrence:")
#print((D2Pnc0/D2Pnc3).simplify().factor())
#print((D2Pnc1/D2Pnc3).simplify().factor())
#print((D2Pnc2/D2Pnc3).simplify().factor())
#print(D2Pnc3)
#print("Second derivative n = 1 recurrence:")
#print((D2P1c0/D2P1c2).simplify().factor())
#print((D2P1c1/D2P1c2).simplify().factor())
#print(D2P1c2.simplify().factor())
#print("Second derivative n = 0 recurrence:")
#print(D2P0c0.simplify().factor())
#print("-------------------------------------------")
## General D3Pn
#D3Pnc0 = (n+a+b)*(n+a+b-1)*(n+a+b-2)*Pnc0/((n+a+b-2)*(n+a+b-3)*(n+a+b-4))
#D3Pnc1 = (n+a+b)*(n+a+b-1)*(n+a+b-2)*Pnc1/((n+a+b-1)*(n+a+b-2)*(n+a+b-3))
#D3Pnc2 = (n+a+b)*(n+a+b-1)*(n+a+b-2)*Pnc2/((n+a+b-1)*(n+a+b-2)*(n+a+b-3))
#D3Pnc3 = 1/(n+a+b-3)
## P1
#D3P1c0 = (a+b+1)*(a+b)*(a+b-1)*P1c0/((a+b)*(a+b-1)*(a+b-2))
#D3P1c1 = (a+b+1)*(a+b)*(a+b-1)*P1c1/((a+b)*(a+b-1)*(a+b-2))
#D3P1c2 = (a+b+1)/(2*(a+b-2))
## P0
#D3P0c0 = 8*(a+b)*(a+b-1)*(a+b-2)
#print("Third derivative recurrence:")
#print((D3Pnc0/D3Pnc3).simplify().factor())
#print((D3Pnc1/D3Pnc3).simplify().factor())
#print((D3Pnc2/D3Pnc3).simplify().factor())
#print(D3Pnc3)
#print("Third derivative n = 1 recurrence:")
#print((D3P1c0/D3P1c2).simplify().factor())
#print((D3P1c1/D3P1c2).simplify().factor())
#print(D3P1c2.simplify().factor())
#print("Third derivative n = 0 recurrence:")
#print(D3P0c0.simplify().factor())
#print("-------------------------------------------")

print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
print("Unity normalization of Worland polynomials")
print("===========================================")
def unitNorm2(n,a,b):
    if n == 0:
        return sy.gamma(a+1)*sy.gamma(b+1)/(2*sy.gamma(a+b+2))
    else:
        return (1/(2*n+a+b+1))*sy.gamma(n+a+1)*sy.gamma(n+b+1)/(2*sy.gamma(n+a+b+1)*sy.gamma(n+1))

# General Pn
UPnc0 = Pnc0*sy.sqrt((unitNorm2(n-2,a,b)/unitNorm2(n,a,b)).simplify())
UPnc1 = Pnc1*sy.sqrt((unitNorm2(n-1,a,b)/unitNorm2(n,a,b)).simplify())
UPnc2 = Pnc2*sy.sqrt((unitNorm2(n-1,a,b)/unitNorm2(n,a,b)).simplify())
UPnc3 = sy.sqrt((2*n+a+b+1)/(n*(n+a)*(n+b)))
# P1
UP1c0 = P1c0*sy.sqrt((unitNorm2(0,a,b)/unitNorm2(1,a,b)).simplify())
UP1c1 = P1c1*sy.sqrt((unitNorm2(0,a,b)/unitNorm2(1,a,b)).simplify())
UP1c2 = sy.Rational(1,2)*sy.sqrt((a+b+3)/((a+1)*(b+1)))
# P0
UP0c0 = P0c0/sy.sqrt(unitNorm2(0,a,b))
print("Polynomial recurrence:")
print(sy.sqrt(((UPnc0/UPnc3)**2).simplify()))
print(sy.sqrt(((UPnc1/UPnc3)**2).simplify()))
print(sy.sqrt(((UPnc2/UPnc3)**2).simplify()))
print(UPnc3)
print("Polynomial n = 1 recurrence:")
print((UP1c0/UP1c2).simplify().factor())
print((UP1c1/UP1c2).simplify().factor())
print((UP1c2).simplify().factor())
print("Polynomial n = 0 recurrence:")
print(UP0c0)
print("-------------------------------------------")

# General first derivative Pn
UDPnc0 = DPnc0*sy.sqrt((unitNorm2(n-1,a-1,b-1)/unitNorm2(n+1,a-1,b-1)).simplify())
UDPnc1 = DPnc1*sy.sqrt((unitNorm2(n,a-1,b-1)/unitNorm2(n+1,a-1,b-1)).simplify())
UDPnc2 = DPnc2*sy.sqrt((unitNorm2(n,a-1,b-1)/unitNorm2(n+1,a-1,b-1)).simplify())
UDPnc3 = sy.sqrt((n+1)*(2*n+a+b+1)/((n+a)*(n+b)))
# DP1
UDP1c0 = DP1c0*sy.sqrt((unitNorm2(1,a-1,b-1)/unitNorm2(2,a-1,b-1)).simplify())
UDP1c1 = DP1c1*sy.sqrt((unitNorm2(1,a-1,b-1)/unitNorm2(2,a-1,b-1)).simplify())
UDP1c2 = (sy.sqrt(2)/2)*sy.sqrt((a+b+1)*(a+b+3)/((a+b)*(a+1)*(b+1)))
# DP0
UDP0c0 = DP0c0/sy.sqrt(unitNorm2(1,a-1,b-1))
print("First derivative recurrence:")
print(sy.sqrt(((UDPnc0/UDPnc3)**2).simplify()))
print(sy.sqrt(((UDPnc1/UDPnc3)**2).simplify()))
print(sy.sqrt(((UDPnc2/UDPnc3)**2).simplify()))
print(UDPnc3)
print("First derivative n = 1 recurrence:")
print(sy.sqrt(((UDP1c0/UDP1c2)**2).simplify()))
print(sy.sqrt(((UDP1c1/UDP1c2)**2).simplify()))
print((UDP1c2).simplify())
print("First derivative n = 0 recurrence:")
print(UDP0c0.simplify())
print("-------------------------------------------")
