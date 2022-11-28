# Black Scholes Merton Model 

Based on the popular European option pricing model [Black-Scholes-Merton.](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model "BSM.")

### Usage
```python
#option = BsmOption(isLong, Type, Spot, Strike, DTE, Risk-Free Rate, Volatility, marketValue, Dividend Yield %)

#Option init using market value
option = BsmOption(True, 'C', 8.86, 10, 18, 0.06, value=2.70)
print("Price = " + str(option.price()))
print("Sigma = " + str(option.sigma))
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
```
