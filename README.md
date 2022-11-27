# Black Scholes Merton Model 

Based on the popular European option pricing model [Black-Scholes-Merton.](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model "BSM.")

### Usage
```python
#option = BsmOption(Spot, Strike, DTE, Risk-Free Rate, Volatility, Dividend Yield %)
option = BsmOption(100, 110, 5, 0.01, 0.4)
print("Price = " + str(option.price('P')))
print("Delta = " + str(option.delta('P')))
print("Gamma = " + str(option.gamma()))
print("Vega  = " + str(option.vega()))
print("Theta = " + str(option.theta('P')))
print("Rho   = " + str(option.rho('P')))
```
