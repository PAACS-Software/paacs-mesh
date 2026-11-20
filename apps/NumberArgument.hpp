// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <charconv>
#include <locale>
#include <sstream>
#include <string>
#include <string_view>

namespace paacs::cli {

template<class T>
bool parseNumber(std::string_view text, T& value) {
    if (text.empty()) return false;

    const char* end = text.data() + text.size();
    if constexpr (requires { std::from_chars(text.data(), end, value); }) {
        const auto [next, error] = std::from_chars(text.data(), end, value);
        return error == std::errc{} && next == end;
    } else {
        // Older libc++ releases provide only the integer from_chars overloads.
        if (text.front() == '+' || text.find_first_of("xX") != std::string_view::npos)
            return false;

        std::istringstream input{std::string(text)};
        input.imbue(std::locale::classic());
        input >> std::noskipws >> value;
        if (!input || !input.eof()) return false;

        // A nonzero decimal that rounded to zero is outside the accepted range.
        const auto significand = text.substr(0, text.find_first_of("eE"));
        return value != T{} || significand.find_first_of("123456789") == std::string_view::npos;
    }
}

} // namespace paacs::cli
