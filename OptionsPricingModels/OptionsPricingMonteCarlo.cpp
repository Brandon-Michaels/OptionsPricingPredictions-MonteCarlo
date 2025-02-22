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

class PolygonAPI {
public:
    static Option fetchOptionData(const std::string& ticker, const std::string& expiration, double strike);
};

size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

Option PolygonAPI::fetchOptionData(const std::string& ticker, const std::string& expiration, double strike) {
    std::string url = "https://api.polygon.io/v3/snapshot/options/" + ticker + "?apiKey=YOUR_API_KEY";
    std::string readBuffer;
    
    CURL* curl = curl_easy_init();
    if (curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        curl_easy_perform(curl);
        curl_easy_cleanup(curl);
    }
    
    Json::Value jsonData;
    Json::Reader reader;
    reader.parse(readBuffer, jsonData);
    
    Option option;
    option.spotPrice = jsonData["underlying_asset"].asDouble();
    option.strikePrice = strike;
    option.timeToMaturity = 0.5; // Assume 6 months until expiration
    option.riskFreeRate = 0.05;
    option.volatility = jsonData["implied_volatility"].asDouble();
    option.type = "call";
    
    return option;
}