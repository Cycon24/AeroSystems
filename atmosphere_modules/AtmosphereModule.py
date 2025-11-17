'''
Atmosphere class that allows for the access of atmoshpheric properties at any height within the atmosophere
'''
import numpy as np
import math
import matplotlib.pyplot as plt

class Atmosphere():
    
    # Selft initialization
    def __init__(self):
        # Create all atmosphereic properties
        Altitude = np.arange(0, 91000,1000) # ft
        
        g = 32.1740484 # ft/s2
        NM = 6076.4   # 1 Nautical Mile = ft
        
        # Define Sea Level values
        T_SL_R = 518.67 # °R
        P_SL = 2116.22  # lbs/ft^2
        rho_SL = 0.00237688  # slugs/ft^3
        a_SL = 1116.45  # ft/s
    
        self.Name = "Plane"
    
        # Atmosheric Array definitions
        self.Alt = Altitude
        self.Temp_R = np.zeros(np.size(Altitude))  # Temperature in Rankin
        self.Temp_F = np.zeros(np.size(Altitude))  # Temperature in Fahrenheit
        
        self.rho = np.zeros(np.size(Altitude))     # Density - slugs/ft^3
        self.P = np.zeros(np.size(Altitude))       # Pressure - lbs/ft^2
        
        self.TR = np.zeros(np.size(Altitude))      # θ Temperature Ratio
        self.PR = np.zeros(np.size(Altitude))      # δ Pressure Ratio
        self.DR = np.zeros(np.size(Altitude))      # σ Density Ratio
        
        self.sqrtDR = np.zeros(np.size(Altitude))  # Square root of Density Ratio
        self.QMS = np.zeros(np.size(Altitude))     # Dynamic Pressure / Mach^2
        self.specW = np.zeros(np.size(Altitude))   # Specific Weight of Air
        self.a = np.zeros(np.size(Altitude))       # Speed of sound ft/s
        self.VELA = np.zeros(np.size(Altitude))    # Speed of sound knots
        self.VRKIN = np.zeros(np.size(Altitude))   #
        
        # Array Calculations
        for i, h in enumerate(Altitude):
            if h <= 36000:
                self.Temp_R[i] = T_SL_R - 3.566*h/1000
                
                self.TR[i] = self.Temp_R[i] / T_SL_R
                self.PR[i] = self.TR[i]**5.2562

            elif h < 82345:
                self.Temp_R[i] = 389.99
                
                self.TR[i] = 0.7519
                self.PR[i] = 0.223361 * math.exp((-0.0481/1000)*(h-36089))
            else:
                self.Temp_R[i] = (- 205.05 + 0.00164*h) + 459.67
                
                self.TR[i] = self.Temp_R[i] / T_SL_R
                self.PR[i] = 0.223361 * math.exp((-0.0481/1000)*(h-36089))
                # self.PR[i] = (51.97/2116) * self.TR[i]**(-11.388)
                
            self.DR[i] = self.PR[i] / self.TR[i]
            self.a[i] = a_SL*np.sqrt(self.TR[i])
            
            self.rho[i] = self.DR[i] * rho_SL
            self.P[i] = self.PR[i] * P_SL
            self.sqrtDR[i] = math.sqrt(self.DR[i])
            self.QMS[i] = 1481.4354 * self.PR[i]
            self.specW[i] = self.rho[i] * g
            self.VELA[i] = 3600 * self.a[i] / NM
            
        self.Temp_F = self.Temp_R - 459.67
        
    
    
    def linterp_h(self, h_range, target_property):
        '''

        Parameters
        ----------
        h_range : Float/Float Array
            Values of altitude in ft to calculate target property at.
        target_property : Atmospheric property array
            The property who's values are needed at the given altitudes.

        Returns
        -------
        Float/Float Array
            A value/array of the target property associated with the altitude inputs.

        '''
        # Linearly interpolates the value of the target property(y) at h(x)
        # Linterp: y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        # Need to get the heights above and below inputted height\
        # can handle vectors
        if type(h_range) == float or type(h_range) == int:
            temp_h = np.empty((1,1), dtype='float')
            temp_h[0] = h_range
            h_range = temp_h
            
            
        propArray = np.zeros(np.size(h_range))
        
        for i_h, h in enumerate(h_range):
            found = False
            for idx, h_test in enumerate(self.Alt):
                if h_test == h:
                    # Values are equal so just return value at h
                    propArray[i_h] = target_property[idx]
                elif h_test >= h:
                    # We have passed the target altitude, so need to break from
                    # loop after saving values
                    h1 = self.Alt[idx-1]
                    h2 = h_test
                    y1 = target_property[idx-1]
                    y2 = target_property[idx]
                    
                    propArray[i_h] = y1 + (y2 - y1) * (h - h1) / (h2 - h1)
                    found = True
                    break
            if not found:
                print('\nError: Target Alitutde not found within Range. At {} ft and {} ft\n'.format(h, h_test))
        
        return propArray if np.size(h_range) > 1 else propArray[0]
    
    @staticmethod
    def linterp(x, x_array, y_array):
        # Linearly interpolates the value of the target property(y) at x
        # Linterp: y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        
        for idx, x_test in enumerate(x_array):
            if x_test == x:
                # Values are equal so just return value at h
                return y_array[idx]
            elif x_test >= x:
                x1 = x_array[idx-1]
                x2 = x_test
                y1 = y_array[idx-1]
                y2 = y_array[idx]
                
                return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        print('\nError: Target Value not found within Range.\n')
        
def StandardAtmosphereTable_Imperial():
    headers = ["h, kft", "P/Pstd", "T/Tstd", "T/Tstd", "T/Tstd", "Z/Tstd", "h, kft"]
    data = """
            0    1.0000  1.0000  0.7708  1.0849  1.0594       0
            1    0.9644  0.9931  0.7972  1.0774  1.0520       1
            2    0.9298  0.9863  0.8237  1.0700  1.0446       2
            3    0.8963  0.9794  0.8501  1.0626  1.0372       3
            4    0.8637  0.9725  0.8575  1.0552  1.0298       4
            5    0.8321  0.9656  0.8575  1.0478  1.0224       5
            6    0.8014  0.9588  0.8575  1.0404  1.0150       6
            7    0.7717  0.9519  0.8575  1.0330  1.0076       7
            8    0.7429  0.9450  0.8575  1.0256  1.0002       8
            9    0.7149  0.9381  0.8575  1.0182  0.9928       9
            10   0.6878  0.9313  0.8565  1.0108  0.9854      10
            11   0.6616  0.9244  0.8502  1.0034  0.9780      11
            12   0.6362  0.9175  0.8438  0.9960  0.9706      12
            13   0.6115  0.9107  0.8375  0.9886  0.9632      13
            14   0.5877  0.9038  0.8312  0.9812  0.9558      14
            15   0.5646  0.8969  0.8248  0.9738  0.9484      15
            16   0.5422  0.8901  0.8185  0.9664  0.9410      16
            17   0.5206  0.8832  0.8121  0.9590  0.9336      17
            18   0.4997  0.8763  0.8058  0.9516  0.9262      18
            19   0.4795  0.8695  0.7994  0.9442  0.9188      19
            20   0.4599  0.8626  0.7931  0.9368  0.9114      20
            21   0.4410  0.8558  0.7867  0.9294  0.9040      21
            22   0.4227  0.8489  0.7804  0.9220  0.8965      22
            23   0.4051  0.8420  0.7740  0.9145  0.8891      23
            24   0.3880  0.8352  0.7677  0.9071  0.8817      24
            25   0.3716  0.8283  0.7613  0.8997  0.8743      25
            """

    # Split the data into lines and then into columns
    lines = data.strip().split('\n')
    data_matrix = [list(map(float, line.split())) for line in lines]
    
    # Convert the list of lists into a NumPy array
    numpy_array = np.array(data_matrix)
    
    # Display the NumPy array
    return(numpy_array)

if __name__=="__main__":
    atm = Atmosphere()
    print(0.85*atm.linterp_h(42975, atm.a))
    print
    # Po / Ps = (1 + gamma - 1)/2 * M^2