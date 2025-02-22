// OptionPricingModels.hpp
#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <curl/curl.h>
#include <json/json.h>
#include <fstream>
#include <string>
#include <cstdlib>

class Option {
public:
    double spotPrice;
    double strikePrice;
    double timeToMaturity;
    double riskFreeRate;
    double volatility;
    std::string type; // "call" or "put"
};

class BlackScholesModel {
public:
    static double calculatePrice(const Option& option);
};

double BlackScholesModel::calculatePrice(const Option& option) {
    double d1 = (log(option.spotPrice / option.strikePrice) +
                 (option.riskFreeRate + 0.5 * option.volatility * option.volatility) * option.timeToMaturity) /
                (option.volatility * sqrt(option.timeToMaturity));
    double d2 = d1 - option.volatility * sqrt(option.timeToMaturity);
    
    double normcdf = [](double x) {
        return 0.5 * erfc(-x * M_SQRT1_2);
    };

    if (option.type == "call") {
        return option.spotPrice * normcdf(d1) - option.strikePrice * exp(-option.riskFreeRate * option.timeToMaturity) * normcdf(d2);
    } else {
        return option.strikePrice * exp(-option.riskFreeRate * option.timeToMaturity) * normcdf(-d2) - option.spotPrice * normcdf(-d1);
    }
}

class MonteCarloModel {
public:
    MonteCarloModel(int numSimulations) : numSimulations_(numSimulations) {}
    double calculatePrice(const Option& option) const;
private:
    int numSimulations_;
};

double MonteCarloModel::calculatePrice(const Option& option) const {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);
    
    double sumPayoff = 0.0;
    double drift = exp((option.riskFreeRate - 0.5 * option.volatility * option.volatility) * option.timeToMaturity);
    double vol = option.volatility * sqrt(option.timeToMaturity);

    #pragma omp parallel for reduction(+:sumPayoff)
    for (int i = 0; i < numSimulations_; ++i) {
        double Z = distribution(generator);
        double St = option.spotPrice * drift * exp(vol * Z);
        double payoff = std::max((option.type == "call" ? St - option.strikePrice : option.strikePrice - St), 0.0);
        sumPayoff += payoff;
    }

    return exp(-option.riskFreeRate * option.timeToMaturity) * (sumPayoff / numSimulations_);
}

class BinomialModel {
public:
    static double calculatePrice(const Option& option);
};

double BinomialModel::calculatePrice(const Option& option) {
    const int steps = 12;
    double dt = option.timeToMaturity / steps;
    double u = exp(option.volatility * sqrt(dt));
    double d = 1 / u;
    double q = (exp(option.riskFreeRate * dt) - d) / (u - d);
    std::vector<double> prices(steps + 1);
    
    for (int i = 0; i <= steps; ++i) {
        prices[i] = std::max((option.type == "call" ? option.spotPrice * pow(u, i) * pow(d, steps - i) - option.strikePrice : option.strikePrice - option.spotPrice * pow(u, i) * pow(d, steps - i)), 0.0);
    }
    
    for (int step = steps - 1; step >= 0; --step) {
        for (int i = 0; i <= step; ++i) {
            prices[i] = exp(-option.riskFreeRate * dt) * (q * prices[i + 1] + (1 - q) * prices[i]);
        }
    }
    
    return prices[0];
}

void loadEnvFile(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Warning: Could not open .env file!" << std::endl;
        return;
    }

    while (std::getline(file, line)) {
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos + 1);
            setenv(key.c_str(), value.c_str(), 1);  // Set environment variable
        }
    }
    file.close();
}

// Option Struct to Store Option Data
struct Option {
    double spotPrice;
    double strikePrice;
    double timeToMaturity;
    double riskFreeRate;
    double volatility;
    std::string type;
};

// Write Callback Function for cURL
size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

// Polygon API Class
class PolygonAPI {
private:
    std::string apiKey;

public:
    // Constructor to Load API Key
    PolygonAPI() {
        loadEnvFile(".env");  // Load API key from .env if available
        const char* key = getenv("POLYGON_API_KEY");

        if (key) {
            apiKey = key;  // Store API key
        } else {
            std::cerr << "Warning: API Key not found! Please set POLYGON_API_KEY in .env or system environment." << std::endl;
            apiKey = "";  // Set empty key (will cause an error if not updated)
        }
    }

    // Fetch Option Data from Polygon.io
    Option fetchOptionData(const std::string& ticker, const std::string& expiration, double strike) {
        if (apiKey.empty()) {
            std::cerr << "Error: API Key is missing. Cannot fetch option data." << std::endl;
            return {};
        }

        std::string url = "https://api.polygon.io/v3/snapshot/options/" + ticker + "?apiKey=" + apiKey;
        std::string readBuffer;

        CURL* curl = curl_easy_init();
        if (curl) {
            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
            curl_easy_perform(curl);
            curl_easy_cleanup(curl);
        } else {
            std::cerr << "Error: Failed to initialize cURL!" << std::endl;
            return {};
        }

        // Parse JSON response
        Json::Value jsonData;
        Json::Reader reader;
        if (!reader.parse(readBuffer, jsonData)) {
            std::cerr << "Error: Failed to parse JSON response!" << std::endl;
            return {};
        }

        // Extract relevant option data
        Option option;
        option.spotPrice = jsonData["underlying_asset"].asDouble();
        option.strikePrice = strike;
        option.timeToMaturity = 0.5; // Assume 6 months to expiration
        option.riskFreeRate = 0.05;
        option.volatility = jsonData["implied_volatility"].asDouble();
        option.type = "call"; // Assuming a call option

        return option;
    }
}
int main() {
    PolygonAPI polygon;

    // Example usage
    std::string ticker = "AAPL";
    std::string expiration = "2024-12-20"; // Expiration date
    double strike = 150.0;

    Option option = polygon.fetchOptionData(ticker, expiration, strike);

    // Display fetched option data
    std::cout << "Spot Price: " << option.spotPrice << std::endl;
    std::cout << "Strike Price: " << option.strikePrice << std::endl;
    std::cout << "Time to Maturity: " << option.timeToMaturity << " years" << std::endl;
    std::cout << "Risk-Free Rate: " << option.riskFreeRate << std::endl;
    std::cout << "Implied Volatility: " << option.volatility << std::endl;
    std::cout << "Option Type: " << option.type << std::endl;

    return 0;
}