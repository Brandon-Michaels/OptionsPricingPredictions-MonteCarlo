import yfinance as yf

# Fetching stock data
def get_stock_data(ticker):
    stock = yf.Ticker(ticker)
    data = stock.history(period="1d")  # Fetch 1 day's data
    return data

# Fetching options data
def get_options_data(ticker, expiration_date):
    stock = yf.Ticker(ticker)
    options = stock.option_chain(expiration_date)  # Fetch option chain for a specific expiration
    calls = options.calls
    puts = options.puts
    return calls, puts

# Example usage
ticker = "AAPL"
stock_data = get_stock_data(ticker)
print(stock_data)

# Get options data for a specific expiration date
expiration_date = "2025-02-28"
calls, puts = get_options_data(ticker, expiration_date)
print("Calls:")
print(calls)
# Fetch AAPL stock data (just an example)
symbol = 'AAPL'
data = yf.download(symbol, period='1d', interval='1m')  # Get minute data for today

# Save the data into a CSV file
data.to_csv('stock_data.csv')

# Option chain example (calls)
option_chain = yf.Ticker(symbol).options
calls = yf.Ticker(symbol).option_chain(option_chain[0]).calls
calls.to_csv('option_data.csv')