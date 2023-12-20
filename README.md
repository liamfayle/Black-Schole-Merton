# Black Scholes Merton Model 

Based on the popular European option pricing model [Black-Scholes-Merton.](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model "BSM.")

### Usage
```python
#option = BsmOption(isLong, Type, Spot, Strike, DTE, Risk-Free Rate, Volatility, marketValue, Dividend Yield %)

#Option init using market value
option = BsmOption(True, 'C', 8.86, 10, 18, 0.06, value=2.70)
print("Price = " + str(option.price()))
print("Sigma = " + str(option.sigma()))
print("Delta = " + str(option.delta()))
print("Gamma = " + str(option.gamma()))
print("Vega  = " + str(option.vega()))
print("Theta = " + str(option.theta()))
print("Rho   = " + str(option.rho()))


#Option init using vol
option = BsmOption(True, 'C', 8.86, 10, 18, 0.06, sigma=.7)
print("Price = " + str(option.price()))
print("Sigma = " + str(option.sigma))
print("Delta = " + str(option.delta()))
print("Gamma = " + str(option.gamma()))
print("Vega  = " + str(option.vega()))
print("Theta = " + str(option.theta()))
print("Rho   = " + str(option.rho()))


#Position init (Short Straddle)
position = OptionPosition()
call = BsmOption(False, 'C', 15.00, 15, 53, 0.05, value=1.75)
put = BsmOption(False, 'P', 15.00, 15, 53, 0.05, value=1.68)
position.addLegs([call, put])
print("Price = " + str(position.price()))
print("Sigma = " + str(position.sigma()))
print("Delta = " + str(position.delta()))
print("Gamma = " + str(position.gamma()))
print("Vega  = " + str(position.vega()))
print("Theta = " + str(position.theta()))
print("Rho   = " + str(position.rho()))
```
