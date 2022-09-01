// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_OPENMP_CONFIG_HPP
#define SEAMAT_OPENMP_CONFIG_HPP

#define SEAMAT_OPENMP_SUPPORT 1

#if defined(SEAMAT_OPENMP_SUPPORT) && (SEAMAT_OPENMP_SUPPORT) == 1
#include <omp.h>
#endif


#endif
