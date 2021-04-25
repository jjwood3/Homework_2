//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "cert-msc51-cpp"
//#pragma ide diagnostic ignored "cert-err58-cpp"
#include <iostream>
#include <random>
#include <cmath>
#include <chrono>

const double STARTING_PRICE = 1868.99;
const double STRIKE_PRICE = 1870.0;
const double VOLATILITY = 0.2979;
const double RISK_FREE_RATE = 0.003866;
const double EXPECTED_DIVIDEND_YIELD = 0.0232;
const double TIME_TO_MATURITY = 1.0/52.0;

const unsigned NUMBER_OF_DIRECT_SIMULATION_RUNS = 6;
const unsigned BASE_NUMBER_OF_REPLICATES_PER_DIRECT_SIMULATION = 1000;

const unsigned NUMBER_OF_ANTITHETIC_SIMULATION_RUNS = 5;
const unsigned BASE_NUMBER_OF_REPLICATES_PER_ANTITHETIC_SIMULATION = 4000;

const unsigned REPLICATE_MULTIPLIER_PER_SIMULATION_RUN = 10;

// #define USE_STATIC_SEED

#ifdef USE_STATIC_SEED
const unsigned STATIC_SEED = 42;
  std::default_random_engine generator( STATIC_SEED );
#else
std::random_device random_device;
std::default_random_engine random_engine( random_device() );
#endif

std::normal_distribution<double> distribution(0.0,1.0);

double generate_next_random_number() {
    return distribution( random_engine );
}

double normal_CDF( double value ) {
    return 0.5 * erfc(-value * M_SQRT1_2);
}

double BSM_call ( double starting_price, double strike_price, double volatility, double risk_free_rate,
                  double time_to_maturity, double expected_dividend_yield ) {

    double d_1 = (
                         log( starting_price / strike_price )
                         + ( risk_free_rate - expected_dividend_yield + 0.5 * pow( volatility, 2 ) ) * time_to_maturity )
                 / ( volatility * sqrt( time_to_maturity )
                 );

    double d_2 = (
                         log( starting_price / strike_price )
                         + ( risk_free_rate - expected_dividend_yield - 0.5 * pow( volatility, 2 ) ) * time_to_maturity )
                 / ( volatility * sqrt( time_to_maturity )
                 );

    return starting_price * exp( -expected_dividend_yield * time_to_maturity ) * normal_CDF( d_1 )
           - strike_price * exp( -risk_free_rate * time_to_maturity ) * normal_CDF( d_2 );
}

void print_bsm_call_price() {
    double bsm_call = BSM_call( STARTING_PRICE, STRIKE_PRICE, VOLATILITY, RISK_FREE_RATE, TIME_TO_MATURITY, EXPECTED_DIVIDEND_YIELD );
    std::cout << "BSM Deterministic Call Price: " << bsm_call << std::endl;
}

void print_csv_header () {
    std::cout << "'sample_size', 'estimated_price', 'estimated_standard_error', '95_percent_confidence_interval_lower',"
                 "'95_percent_confidence_interval_upper, 'runtime_in_seconds', 'efficiency'" << std::endl;
}

void print_simulation_results( double replicate_count, double estimated_price, double estimated_standard_error,
                               double confidence_interval_lower, double confidence_interval_upper, double runtime_in_seconds ) {
    std::cout << replicate_count << ", ";
    std::cout << estimated_price << ", ";
    std::cout << estimated_standard_error << ", ";
    std::cout << confidence_interval_lower << ", ";
    std::cout << confidence_interval_upper << ", ";
    std::cout << runtime_in_seconds << ", ";
    std::cout << estimated_standard_error * estimated_standard_error * runtime_in_seconds << std::endl;
}

const double BSM_DETERMINISTIC_PART = STARTING_PRICE * exp( ( RISK_FREE_RATE - EXPECTED_DIVIDEND_YIELD - 0.5 * VOLATILITY * VOLATILITY ) * TIME_TO_MATURITY );
const double BSM_RANDOM_COEFFICIENT = VOLATILITY * sqrt( TIME_TO_MATURITY );

double bma_deterministic_part = STARTING_PRICE * exp(
        ( RISK_FREE_RATE - EXPECTED_DIVIDEND_YIELD - 0.5 * VOLATILITY * VOLATILITY ) * TIME_TO_MATURITY
);

double bma_stocastic_part_coefficient = VOLATILITY * sqrt( TIME_TO_MATURITY );

double simulate_new_bma_price_direct() {
    double Si = bma_deterministic_part * exp( bma_stocastic_part_coefficient * generate_next_random_number() );
    return ( Si - STRIKE_PRICE > 0 )
           ? exp( -1 * RISK_FREE_RATE * TIME_TO_MATURITY ) * ( Si - STRIKE_PRICE )
           : 0;
}

double simulate_new_bma_price_antithetic () {
    double random_number = generate_next_random_number();
    double Si_1 = bma_deterministic_part * exp( bma_stocastic_part_coefficient * random_number );
    double Si_2 = bma_deterministic_part * exp( -1.0 * bma_stocastic_part_coefficient * random_number );
    Si_1 = ( Si_1 - STRIKE_PRICE > 0 )
           ? ( Si_1 - STRIKE_PRICE )
           : 0;
    Si_2 = ( Si_2 - STRIKE_PRICE > 0 )
           ? ( Si_2 - STRIKE_PRICE )
           : 0;
    return 0.5 * exp( -1 * RISK_FREE_RATE * TIME_TO_MATURITY ) * ( Si_1 + Si_2 );
}

void run_direct_simulation( int simulation_number ) {
    auto start = std::chrono::high_resolution_clock::now();
    int replicate_count = BASE_NUMBER_OF_REPLICATES_PER_DIRECT_SIMULATION * pow( REPLICATE_MULTIPLIER_PER_SIMULATION_RUN, simulation_number );
    double x_1 = simulate_new_bma_price_direct();
    double x_bar = x_1;
    double y_bar = x_1 * x_1;
    for( int k = 2; k <= replicate_count; k++ ) {
        double x_k = simulate_new_bma_price_direct();
        x_bar = ( 1 - 1.0/k ) * x_bar + 1.0/k * x_k;
        y_bar = ( 1 - 1.0/k ) * y_bar + 1.0/k * x_k * x_k;
    }
    double standard_error = sqrt( 1.0 / ( replicate_count ) * ( y_bar - x_bar * x_bar ) );
    double confidence_interval_lower  = x_bar - 1.96 * VOLATILITY / sqrt( replicate_count );
    double confidence_interval_upper  = x_bar + 1.96 * VOLATILITY / sqrt( replicate_count );
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    double elapsed_time = elapsed.count();
    print_simulation_results( replicate_count, x_bar, standard_error, confidence_interval_lower, confidence_interval_upper, elapsed_time );
}

void run_antithetic_simulation( int simulation_number ) {
    auto start = std::chrono::high_resolution_clock::now();
    int replicate_count = BASE_NUMBER_OF_REPLICATES_PER_ANTITHETIC_SIMULATION * pow( REPLICATE_MULTIPLIER_PER_SIMULATION_RUN, simulation_number );
    double x_1 = simulate_new_bma_price_antithetic();
    double x_bar = x_1;
    double y_bar = x_1 * x_1;
    for( int k = 2; k <= replicate_count; k++ ) {
        double x_k = simulate_new_bma_price_antithetic();
        x_bar = ( 1 - 1.0/k ) * x_bar + 1.0/k * x_k;
        y_bar = ( 1 - 1.0/k ) * y_bar + 1.0/k * x_k * x_k;
    }
    double standard_error = sqrt( 1.0 / ( replicate_count ) * ( y_bar - x_bar * x_bar ) );
    double confidence_interval_lower  = x_bar - 1.96 * VOLATILITY / sqrt( replicate_count );
    double confidence_interval_upper  = x_bar + 1.96 * VOLATILITY / sqrt( replicate_count );
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    double elapsed_time = elapsed.count();
    print_simulation_results( replicate_count, x_bar, standard_error, confidence_interval_lower, confidence_interval_upper, elapsed_time );
}

int main() {
    print_bsm_call_price();

    std::cout << std::endl;
    std::cout << "CSV Data table for Stochastic SBM Simulation using Direct Method." << std::endl;
    print_csv_header();
    for( int i = 0; i < NUMBER_OF_DIRECT_SIMULATION_RUNS; i++ ) {
        run_direct_simulation(i);
    }

    std::cout << std::endl;
    std::cout << "CSV Data table for Stochastic SBM Simulation using Antithetic Method." << std::endl;
    print_csv_header();
    for( int i = 0; i < NUMBER_OF_ANTITHETIC_SIMULATION_RUNS; i++ ) {
        run_antithetic_simulation(i);
    }
    return 0;
}

//#pragma clang diagnostic pop
