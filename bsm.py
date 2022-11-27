from scipy.stats import norm
import numpy as np


class BsmOption:
    def __init__(self, S, K, T, r, sigma, q=0):
        '''
        S -> Underlying Price           [$]             [100]           100$ Underlying                 \n
        K -> Strike                     [$]             [110]           110$ Strike                     \n
        T -> Time until expiration      [Decimal]       [20]            20 DTE                          \n
        r -> Risk free rate             [Decimal]       [0.01]          1% RFR Continous yield          \n
        sigma -> Volatility             [Decimal]       [0.45]          45% Vol                         \n
        q -> Dividend Value             [Decimal]       [0.01]          1% Continous Div yield            
        '''
        self.S = S 
        self.K = K
        self.T = T / 365
        self.r = r
        self.sigma = sigma
        self.q = q

    
    @staticmethod
    def N(x):
        return norm.cdf(x)

    @staticmethod
    def N_prime(x):
        return norm.pdf(x)

    @property
    def params(self):
        return {'S': self.S,
                'K': self.K,
                'T': self.T,
                'r': self.r,
                'sigma': self.sigma,
                'q': self.q}


    def d1(self):
        return (np.log(self.S/self.K) + ((self.r -self.q + self.sigma**2/2)*self.T)) / (self.sigma*np.sqrt(self.T))
    
    def d2(self):
        return self.d1() - self.sigma*np.sqrt(self.T)
    
    def _call_value(self):
        return self.S*np.exp(-self.q*self.T)*self.N(self.d1()) - self.K*np.exp(-self.r*self.T) * self.N(self.d2())
                    
    def _put_value(self):
        return self.K*np.exp(-self.r*self.T) * self.N(-self.d2()) - self.S*np.exp(-self.q*self.T)*self.N(-self.d1())

    def delta(self, type_):
        '''
        Return Delta Greek Value \n
        C -> CAll \n
        P -> PUT \n
        B -> BOTH \n
        '''
        if type_ == 'C':
            return self.N(self.d1())
        if type_ == 'P':
            return self.delta('C') - 1
        if type_ == 'B':
            return  {'call': self.delta('C'), 'put': self.delta('P')}
        else:
            raise ValueError('Unrecognized type')
        

    def gamma(self):
        '''
        Return Gamma Greek Value \n
        '''
        return (self.N_prime(self.d1())) / (self.S * self.sigma * np.sqrt(self.T))

    def vega(self):
        '''
        Return Delta Greek Value \n
        '''
        return self.S * self.N_prime(self.d1()) * np.sqrt(self.T)


    def theta(self, type_):
        '''
        Return Delta Greek Value \n
        C -> CAll \n
        P -> PUT \n
        B -> BOTH \n
        '''
        if type_ == 'C':
            return ( - (self.S * self.N_prime(self.d1()) * self.sigma) / (2 * np.sqrt(self.T)) ) - (self.r * self.K * np.exp(-self.r * self.T) * self.N(self.d2()))
        if type_ == 'P':
            return ( - (self.S * self.N_prime(self.d1()) * self.sigma) / (2 * np.sqrt(self.T)) ) + (self.r * self.K * np.exp(-self.r * self.T) * self.N(-self.d2()))
        if type_ == 'B':
            return  {'call': self.theta('C'), 'put': self.theta('P')}
        else:
            raise ValueError('Unrecognized type')


    def rho(self, type_):
        '''
        Return rho Greek Value \n
        C -> CAll \n
        P -> PUT \n
        B -> BOTH \n
        '''
        if type_ == 'C':
            return self.K * self.T * np.exp(-self.r * self.T) * self.N(self.d2())
        if type_ == 'P':
            return -self.K * self.T * np.exp(-self.r * self.T) * self.N(-self.d2())
        if type_ == 'B':
            return  {'call': self.rho('C'), 'put': self.rho('P')}
        else:
            raise ValueError('Unrecognized type')
    
    def price(self, type_ = 'C'):
        '''
        C -> CAll \n
        P -> PUT \n
        B -> BOTH \n
        '''
        if type_ == 'C':
            return self._call_value()
        if type_ == 'P':
            return self._put_value() 
        if type_ == 'B':
            return  {'call': self._call_value(), 'put': self._put_value()}
        else:
            raise ValueError('Unrecognized type')


    def set_sigma(self, sigma):
        self.sigma = sigma


