#include "../common/consts.hpp"
#include <iostream>

template<typename F>
unsigned find_count_pade_members() {
    F ln2_n = 1;
    for (unsigned n = 1; n < 1000; n++) {
        ln2_n *= adaai::C_LN_2<F>;
        if (ln2_n < adaai::C_EPS<F>) {
            return n - 1;
        }
    }
    return -1;
}

int main() {
    /*
     Знаем что ошибка в аппроксимация Паде = O(x^{n+m+1})
     т.к. |x| < ln2, то для каждого типа можем посчитать, сколько членов многочлена Паде нужно взять
     */
    std::cout << "float value n + m = " << find_count_pade_members<float>() << "\n";
    std::cout << "double value n + m = " << find_count_pade_members<double>() << "\n";
    std::cout << "long double value n + m = " << find_count_pade_members<long double>() << "\n";

}