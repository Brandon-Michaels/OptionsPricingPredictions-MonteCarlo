#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <random>

// CDF of standard normal distribution
double normCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// Black-Scholes Model
double blackScholes(double S, double K, double T, double r, double sigma, bool callOption) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);

    double callPrice = S * normCDF(d1) - K * exp(-r * T) * normCDF(d2);
    double putPrice = K * exp(-r * T) * normCDF(-d2) - S * normCDF(-d1);

    return callOption ? callPrice : putPrice;
}

// Binomial Tree Model
double binomialOptionPricing(double S, double K, double T, double r, double sigma, int N, bool callOption) {
    double dt = T / N;
    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p = (exp(r * dt) - d) / (u - d);
    std::vector<double> prices(N + 1, 0.0);

    // Calculate option prices at final nodes
    for (int i = 0; i <= N; i++) {
        double ST = S * pow(u, N - i) * pow(d, i);
        prices[i] = std::max(callOption ? (ST - K) : (K - ST), 0.0);
    }

    // Step backward to calculate the option price at the root node
    for (int j = N - 1; j >= 0; j--) {
        for (int i = 0; i <= j; i++) {
            prices[i] = exp(-r * dt) * (p * prices[i] + (1 - p) * prices[i + 1]);
        }
    }

    return prices[0];
}

// Monte Carlo Simulation Model
double monteCarloOptionPricing(double S, double K, double T, double r, double sigma, int numSimulations, bool callOption) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);

    double totalPayoff = 0.0;
    for (int i = 0; i < numSimulations; i++) {
        double ST = S;
        for (int j = 0; j < T * 252; j++) {  // Assume daily trading (252 trading days in a year)
            ST *= exp((r - 0.5 * sigma * sigma) / 252 + sigma * d(gen) * sqrt(1.0 / 252));
        }
        totalPayoff += std::max(callOption ? (ST - K) : (K - ST), 0.0);
    }

    return exp(-r * T) * totalPayoff / numSimulations;
}

// Helper function to parse the CSV for stock data
void readCSV(std::vector<std::vector<std::string> >& data, const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }
        data.push_back(row);
    }
}

int main() {
    // Read data from CSV (stock data)
    std::vector<std::vector<std::string> > stockData;
    readCSV(stockData, "stock_data.csv");

    // Read option chain data (calls)
    std::vector<std::vector<std::string> > optionData;
    readCSV(optionData, "option_data.csv");

    if (stockData.empty() || optionData.empty()) {
        std::cerr << "Error: CSV data is empty." << std::endl;
        return 1;
    }

    try {
        // Ensure we have enough rows and columns before accessing
        if (stockData.size() < 2 || stockData.back().size() <= 4) {
            throw std::invalid_argument("Stock data missing or incorrectly formatted.");
        }
        if (optionData.size() < 2 || optionData[0].size() <= 6) {
            throw std::invalid_argument("Option data missing or incorrectly formatted.");
        }

        // Convert necessary values with additional checks
        double stockPrice = std::stod(stockData.back()[4]);  // Latest stock closing price
        double strikePrice = std::stod(optionData[1][2]);    // First option strike price
        double sigma = std::stod(optionData[1][6]);          // Implied volatility

        double T = 30.0 / 365;  // 30 days to expiration
        double r = 0.05;        // Risk-free rate
        bool callOption = true; // Call option

        // Black-Scholes model
        double bsPrice = blackScholes(stockPrice, strikePrice, T, r, sigma, callOption);
        std::cout << "Black-Scholes Price: " << bsPrice << std::endl;

        // Binomial model
        double binomialPrice = binomialOptionPricing(stockPrice, strikePrice, T, r, sigma, 1000, callOption);
        std::cout << "Binomial Model Price: " << binomialPrice << std::endl;

        // Monte Carlo simulation
        double monteCarloPrice = monteCarloOptionPricing(stockPrice, strikePrice, T, r, sigma, 10000, callOption);
        std::cout << "Monte Carlo Simulation Price: " << monteCarloPrice << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
