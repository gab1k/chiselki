#ifndef RANDOM_NUM_HPP_
#define RANDOM_NUM_HPP_

#include <random>

template<typename T>
T get_rand_val(T from, T to) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<T> distribution(from, to);
    return distribution(generator);
}

#endif // RANDOM_NUM_HPP_