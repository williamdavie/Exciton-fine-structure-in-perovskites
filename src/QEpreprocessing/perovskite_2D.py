'''
Base class for 2D pervoskite structure, general formula (A)2(BX4)
'''

class perovskite2D:

    
    def __init__(self, B: str, Bmass: float, X: str, Xmass: float, A: str='Cs', Amass: float=132.9):
        
        self.B = B
        self.Bmass = Bmass
        self.X = X
        self.Xmass = Xmass
        self.A = A
        self.Amass = Amass