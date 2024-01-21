# simulator.py is a program that simulates the combustion of a solid rocket motor
# The program calculates the thrust, chamber pressure, and other parameters of the motor
# The program also plots the thrust, chamber pressure, and other parameters of the motor
# Writen by: Manfred Gawlas PWRINSPACE

import math
import matplotlib.pyplot as plt

# Class to store cylinder geometry and burn distance
class AbGeometry:
    def __init__(self, L, R, D, x):
        self.L = L
        # Radius port
        self.R = R
        self.D = D
        # Burn distance
        self.x = x

# Class to store general parameters
class General:
    def __init__(self, Pch, L, D, d, At, Ae, Pe):
        self.Pch = Pch
        self.L = L
        self.D = D
        self.d = d
        self.At = At
        self.Ae = Ae
        self.Pe = Pe

# Class to store fuel parameters
class Fuel:
    def __init__(self, Cstar, a, n, k, rho):
        self.Cstar = Cstar
        self.a = a
        self.n = n
        self.k = k
        # Propellant density
        self.rho = rho

# Function to calculate CF
def CF(Cylinder, fuel):
    return math.pow((2*fuel.k*fuel.k)/(fuel.k-1) * math.pow(2/(fuel.k+1), (fuel.k+1)/(fuel.k-1)) * (1 - math.pow(Cylinder.Pe/Cylinder.Pch, (fuel.k-1)/fuel.k)), 0.5)

# Function to calculate Thrust
def F(Cylinder, fuel):
    return CF(Cylinder, fuel) * Cylinder.At * Cylinder.Pch + (Cylinder.Pe - 101325) * Cylinder.Ae

# Function to calculate Area burning
def FunctionAb(Ab):
    return 2 * 3.14 * ((Ab.R + Ab.x) * (Ab.L - 2*Ab.x) + (Ab.D/2) * (Ab.D/2) - (Ab.R + Ab.x) * (Ab.R + Ab.x))


#############################################################################################################
# Different method to calculate Pe (not used)
# Idea is to calculate Pe for different values of Pch and store them in a list, and then use the list to find Pe
# This method is not used because for certain values of Pch, Pe is not calculated correctly (Pe is too high) and is less time efficient
# However, this method is way cooler as an idea
#############################################################################################################
'''
def compute_Pe(Cylinder, fuel):
    Pe_tab = []
    Pe_tab.append(0)
    for i in range(1, 5000):
        output_num = num_Pe(Cylinder, fuel, float(i))
        if(output_num > 3 * Pe_tab[i-1] and i>10):
            Pe_tab.append(Pe_tab[i-1])
        else:
            Pe_tab.append(output_num)
    return Pe_tab
    
def num_Pe(Cylinder, fuel, Pch):
    Pe=float(0)
    Pch = float(Pch * 1000)
    min = 100000
    Pe_out = 0
    E = Cylinder.At/Cylinder.Ae     
    while(Pe<Pch):
        Pe+=1000
        temp = abs(math.pow((fuel.k +1)/2, 1/(fuel.k-1)) * math.pow(Pe/Pch, 1/fuel.k) * math.sqrt((fuel.k+1)/(fuel.k-1) * (1 - math.pow(Pe/Pch, 1/(fuel.k-1)))) - E)
        
        if(temp<min):
            min=temp
            Pe_out=Pe
    return Pe_out
'''    
#############################################################################################################


# Function to calculate Pe numerically
def num_comp_Pe(Cylinder, fuel):
    Pe=float(101325-100000) # Initial value of Pe, it's not 0, to avoid uselles calculations
    min = 1000 # Minimal value to compare to
    Pe_out = float(0) # Output value of Pe
    E = Cylinder.At/Cylinder.Ae # (Expansion ratio)^(-1)

    # While loop to calculate Pe
    while(Pe<101325+50000): # Works thill certaion value of Pe, that I choose since higher values are not likely to happen
        Pe+=100 # Increment Pe

        # Calculate delta E = E(Pe) - E(const), where E(const) is E inputed and E(Pe) is calculated for the current Pe
        temp = abs(math.pow((fuel.k +1)/2, 1/(fuel.k-1)) * math.pow(Pe/Cylinder.Pch, 1/fuel.k) * math.sqrt((fuel.k+1)/(fuel.k-1) * (1 - math.pow(Pe/Cylinder.Pch, (fuel.k-1)/fuel.k))) - E)
        
        # If delta E is smaller than the minimal value, then delta E is the new minimal value and Pe is the new output value
        if(temp<min):
            min=temp
            Pe_out=Pe

    return Pe_out

# Main function to calculate all data of simulation
def calculate_parameters():
    # Initial values and definitions
    r=0
    Ic=0
    F0=0
    tc=0
    Deltat=0.001

    # Object definitions
    Cylinder = General(0, 0.12, 0.051, 0.03, (0.01/2)*(0.01/2)*3.14, (0.014/2)*(0.014/2)*3.14, 101325)
    fuel = Fuel(779, 0.000016, 0.371, 1.18, 1848)
    Ab = AbGeometry(Cylinder.L, Cylinder.d/2, Cylinder.D, 0)

    # Calculate initial value of Pch
    Cylinder.Pch=math.pow(FunctionAb(Ab) / Cylinder.At, 1/(1-fuel.n)) * math.pow(fuel.Cstar * fuel.rho * fuel.a, 1/(1-fuel.n))

    #Pe_tab = compute_Pe(Cylinder, fuel)

    F0_list = []  # List to store F0 values
    tc_list = []  # List to store tc values
    Pch_list = []  # List to store Pch values
    Kn_list = []  # List to store Kn values
    Pe_list = []  # List to store Pe values

    # While loop to calculate simulation data
    while(Ab.x<((Cylinder.D-Cylinder.d)/2)):
        #Cylinder.Pe = Pe_tab[int(Cylinder.Pch/1000)]

        # Calculate Pe numerically
        Cylinder.Pe = num_comp_Pe(Cylinder, fuel)

        
        tc+=Deltat # Increment total time
        F0=F(Cylinder, fuel) # Calculate thrust
        Ic+=F0 * Deltat # Increment impulse

        r=fuel.a * math.pow(Cylinder.Pch, fuel.n) # Calculate regression rate
        Ab.x+=r * Deltat # Increment burn distance

        # Calculate new Pch
        Cylinder.Pch=math.pow(FunctionAb(Ab) / Cylinder.At, 1/(1-fuel.n)) * math.pow(fuel.Cstar * fuel.rho * fuel.a, 1/(1-fuel.n))

        F0_list.append(F0)  # Append F0 value to the list
        tc_list.append(tc)  # Append tc value to the list
        Pch_list.append(Cylinder.Pch)  # Append Pch value to the list
        Kn_list.append(FunctionAb(Ab) / Cylinder.At)  # Append Kn value to the list
        Pe_list.append(Cylinder.Pe)  # Append Pe value to the list

    return Ic, tc, Cylinder.Pch, Cylinder.At, F0_list, tc_list, Pch_list, Kn_list, Pe_list

# Function to print data for the user
def ShowOutput():
    Ic, tc, Pch, At, F0_list, tc_list, Pch_list, Kn_list, Pe_list = calculate_parameters()

    print("Ic=", Ic)
    print("tc=", tc)
    print("Pch=", Pch)
    print("At=", At)

    # Plot F0 vs Time
    plt.plot(tc_list, F0_list)
    plt.xlabel('Time [s]')
    plt.ylabel('Thrust [N]')
    plt.title('Thrust vs Time')
    plt.show()

    # Plot Pch vs Time
    plt.plot(tc_list, Pch_list)
    plt.xlabel('Time [s]')
    plt.ylabel('Chamber Pressure [Pa]')
    plt.title('Chamber Pressure vs Time')
    plt.show()

    # Plot Kn vs Time
    plt.plot(tc_list, Kn_list)
    plt.xlabel('Time [s]')
    plt.ylabel('Kn [ ]')
    plt.title('Kn vs Time')
    plt.show()

    # Plot Pch vs Pe
    plt.plot(Pch_list, Pe_list)
    plt.xlabel('Chamber Pressure [Pa]')
    plt.ylabel('Exit Pressure [Pa]')
    plt.title('Exit Pressure vs Chamber Pressure')
    plt.show()

ShowOutput()