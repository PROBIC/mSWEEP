// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_MATH_UTIL_HPP
#define SEAMAT_MATH_UTIL_HPP

#include <cmath>
#include <limits>

namespace seamat {
    template <typename T>
    bool nearly_equal(const T a, const T b) {
	return std::nextafter(a, std::numeric_limits<T>::lowest()) <= b
	       && std::nextafter(a, std::numeric_limits<T>::max()) >= b;
    }

    template <typename T>
    bool nearly_equal(const T a, const T b, const uint64_t factor) {
	/* a factor of epsilon */
	T min_a = a - (a - std::nextafter(a, std::numeric_limits<T>::lowest())) * factor;
	T max_a = a + (std::nextafter(a, std::numeric_limits<T>::max()) - a) * factor;

	return min_a <= b && max_a >= b;
    }
}

#endif
