#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

// Do not use "using namespace std;"

namespace sjtu {
class int2048 {
private:
  std::vector<int> digits; // stored in reverse order (least significant first)
  bool negative;
  
  void remove_leading_zeros();
  static int compare_abs(const int2048 &a, const int2048 &b);
  static int2048 add_abs(const int2048 &a, const int2048 &b);
  static int2048 sub_abs(const int2048 &a, const int2048 &b);
  void divide_by_2();
  static int2048 multiply_fft(const int2048 &a, const int2048 &b);
  static void fft(std::vector<std::complex<double>> &a, bool invert);
  
public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};

// Base for digit compression - using 10^9 for efficiency
static const long long BASE = 1000000000LL;
static const int BASE_DIGITS = 9;

// Private helper methods
void int2048::remove_leading_zeros() {
    while (digits.size() > 1 && digits.back() == 0) {
        digits.pop_back();
    }
    if (digits.size() == 1 && digits[0] == 0) {
        negative = false;
    }
}

int int2048::compare_abs(const int2048 &a, const int2048 &b) {
    if (a.digits.size() != b.digits.size()) {
        return a.digits.size() < b.digits.size() ? -1 : 1;
    }
    for (int i = a.digits.size() - 1; i >= 0; --i) {
        if (a.digits[i] != b.digits[i]) {
            return a.digits[i] < b.digits[i] ? -1 : 1;
        }
    }
    return 0;
}

int2048 int2048::add_abs(const int2048 &a, const int2048 &b) {
    int2048 result;
    result.negative = false;
    int carry = 0;
    size_t max_size = std::max(a.digits.size(), b.digits.size());
    result.digits.resize(max_size);
    
    for (size_t i = 0; i < max_size || carry; ++i) {
        if (i == result.digits.size()) {
            result.digits.push_back(0);
        }
        long long sum = carry;
        if (i < a.digits.size()) sum += a.digits[i];
        if (i < b.digits.size()) sum += b.digits[i];
        result.digits[i] = sum % BASE;
        carry = sum / BASE;
    }
    return result;
}

int2048 int2048::sub_abs(const int2048 &a, const int2048 &b) {
    int2048 result;
    result.negative = false;
    result.digits = a.digits;
    
    int carry = 0;
    for (size_t i = 0; i < b.digits.size() || carry; ++i) {
        long long diff = result.digits[i] - carry;
        if (i < b.digits.size()) {
            diff -= b.digits[i];
        }
        carry = 0;
        if (diff < 0) {
            diff += BASE;
            carry = 1;
        }
        result.digits[i] = diff;
    }
    result.remove_leading_zeros();
    return result;
}

void int2048::divide_by_2() {
    int carry = 0;
    for (int i = digits.size() - 1; i >= 0; --i) {
        long long cur = (long long)carry * BASE + digits[i];
        digits[i] = cur / 2;
        carry = cur % 2;
    }
    remove_leading_zeros();
}

void int2048::fft(std::vector<std::complex<double>> &a, bool invert) {
    int n = a.size();
    if (n == 1) return;
    
    // Bit reversal
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * M_PI / len * (invert ? -1 : 1);
        std::complex<double> wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (int j = 0; j < len / 2; ++j) {
                std::complex<double> u = a[i + j];
                std::complex<double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    
    if (invert) {
        for (auto &x : a) {
            x /= n;
        }
    }
}

int2048 int2048::multiply_fft(const int2048 &a, const int2048 &b) {
    typedef std::complex<double> base;
    
    int result_size = 1;
    while (result_size < (int)(a.digits.size() + b.digits.size())) {
        result_size <<= 1;
    }
    
    std::vector<base> fa(a.digits.begin(), a.digits.end());
    std::vector<base> fb(b.digits.begin(), b.digits.end());
    fa.resize(result_size);
    fb.resize(result_size);
    
    fft(fa, false);
    fft(fb, false);
    
    for (int i = 0; i < result_size; ++i) {
        fa[i] *= fb[i];
    }
    
    fft(fa, true);
    
    int2048 result;
    result.digits.resize(result_size);
    long long carry = 0;
    for (int i = 0; i < result_size; ++i) {
        long long cur = (long long)(fa[i].real() + 0.5) + carry;
        result.digits[i] = cur % BASE;
        carry = cur / BASE;
    }
    while (carry) {
        result.digits.push_back(carry % BASE);
        carry /= BASE;
    }
    
    result.remove_leading_zeros();
    return result;
}

// Constructors
int2048::int2048() : negative(false) {
    digits.push_back(0);
}

int2048::int2048(long long x) {
    if (x < 0) {
        negative = true;
        x = -x;
    } else {
        negative = false;
    }
    if (x == 0) {
        digits.push_back(0);
    } else {
        while (x > 0) {
            digits.push_back(x % BASE);
            x /= BASE;
        }
    }
}

int2048::int2048(const std::string &s) {
    read(s);
}

int2048::int2048(const int2048 &other) : digits(other.digits), negative(other.negative) {}

// Read from string
void int2048::read(const std::string &s) {
    digits.clear();
    negative = false;
    
    int start = 0;
    if (s[0] == '-') {
        negative = true;
        start = 1;
    } else if (s[0] == '+') {
        start = 1;
    }
    
    // Parse from right to left in chunks of BASE_DIGITS
    for (int i = s.length(); i > start; i -= BASE_DIGITS) {
        int left = std::max(start, i - BASE_DIGITS);
        int digit = 0;
        for (int j = left; j < i; ++j) {
            digit = digit * 10 + (s[j] - '0');
        }
        digits.push_back(digit);
    }
    
    remove_leading_zeros();
}

// Print
void int2048::print() {
    if (negative) std::cout << '-';
    std::cout << digits.back();
    for (int i = digits.size() - 2; i >= 0; --i) {
        // Pad with leading zeros
        int d = digits[i];
        for (int j = BASE_DIGITS - 1; j > 0; --j) {
            if (d < (int)pow(10, j)) std::cout << '0';
            else break;
        }
        std::cout << d;
    }
}

// Add operations
int2048 &int2048::add(const int2048 &other) {
    bool original_sign = negative;
    if (negative == other.negative) {
        *this = add_abs(*this, other);
        this->negative = original_sign;
    } else {
        int cmp = compare_abs(*this, other);
        if (cmp >= 0) {
            *this = sub_abs(*this, other);
            this->negative = original_sign;
        } else {
            *this = sub_abs(other, *this);
            this->negative = other.negative;
        }
    }
    return *this;
}

int2048 add(int2048 a, const int2048 &b) {
    a.add(b);
    return a;
}

// Subtract operations
int2048 &int2048::minus(const int2048 &other) {
    bool original_sign = negative;
    if (negative != other.negative) {
        *this = add_abs(*this, other);
        this->negative = original_sign;
    } else {
        int cmp = compare_abs(*this, other);
        if (cmp == 0) {
            *this = int2048(0);
        } else if (cmp > 0) {
            *this = sub_abs(*this, other);
            this->negative = original_sign;
        } else {
            *this = sub_abs(other, *this);
            this->negative = !original_sign;
        }
    }
    return *this;
}

int2048 minus(int2048 a, const int2048 &b) {
    a.minus(b);
    return a;
}

// Unary operators
int2048 int2048::operator+() const {
    return *this;
}

int2048 int2048::operator-() const {
    int2048 result = *this;
    if (!(result.digits.size() == 1 && result.digits[0] == 0)) {
        result.negative = !result.negative;
    }
    return result;
}

// Assignment
int2048 &int2048::operator=(const int2048 &other) {
    if (this != &other) {
        digits = other.digits;
        negative = other.negative;
    }
    return *this;
}

// Addition
int2048 &int2048::operator+=(const int2048 &other) {
    return add(other);
}

int2048 operator+(int2048 a, const int2048 &b) {
    a += b;
    return a;
}

// Subtraction
int2048 &int2048::operator-=(const int2048 &other) {
    return minus(other);
}

int2048 operator-(int2048 a, const int2048 &b) {
    a -= b;
    return a;
}

// Multiplication
int2048 &int2048::operator*=(const int2048 &other) {
    bool result_negative = negative != other.negative;
    
    // Simple multiplication for small numbers
    if (digits.size() + other.digits.size() < 100) {
        std::vector<long long> result(digits.size() + other.digits.size(), 0);
        for (size_t i = 0; i < digits.size(); ++i) {
            for (size_t j = 0; j < other.digits.size(); ++j) {
                result[i + j] += (long long)digits[i] * other.digits[j];
            }
        }
        
        digits.clear();
        long long carry = 0;
        for (size_t i = 0; i < result.size() || carry; ++i) {
            if (i >= result.size()) result.push_back(0);
            result[i] += carry;
            carry = result[i] / BASE;
            digits.push_back(result[i] % BASE);
        }
    } else {
        // For very large numbers, use FFT-based multiplication
        *this = multiply_fft(*this, other);
    }
    
    negative = result_negative;
    remove_leading_zeros();
    return *this;
}

int2048 operator*(int2048 a, const int2048 &b) {
    a *= b;
    return a;
}

// Division using long division algorithm
int2048 &int2048::operator/=(const int2048 &other) {
    bool result_negative = negative != other.negative;
    
    int2048 dividend = *this;
    dividend.negative = false;
    int2048 divisor = other;
    divisor.negative = false;
    
    int cmp = compare_abs(dividend, divisor);
    if (cmp < 0) {
        // Result is 0 or -1 depending on signs
        if (result_negative) {
            *this = int2048(-1);
        } else {
            *this = int2048(0);
        }
        return *this;
    }
    
    if (cmp == 0) {
        *this = int2048(result_negative ? -1 : 1);
        return *this;
    }
    
    // Long division algorithm
    int2048 quotient(0);
    int2048 remainder(0);
    
    // Process from most significant to least significant digit
    for (int i = dividend.digits.size() - 1; i >= 0; --i) {
        // Shift remainder left by one base position and add current digit
        remainder.digits.insert(remainder.digits.begin(), dividend.digits[i]);
        remainder.remove_leading_zeros();
        
        if (compare_abs(remainder, divisor) < 0) {
            quotient.digits.insert(quotient.digits.begin(), 0);
            continue;
        }
        
        // Binary search for the quotient digit
        int left = 0, right = BASE - 1;
        while (left < right) {
            int mid = (left + right + 1) / 2;
            int2048 temp = divisor;
            // Multiply divisor by mid
            long long carry = 0;
            for (size_t j = 0; j < temp.digits.size(); ++j) {
                long long prod = (long long)temp.digits[j] * mid + carry;
                temp.digits[j] = prod % BASE;
                carry = prod / BASE;
            }
            if (carry > 0) {
                temp.digits.push_back(carry);
            }
            
            if (compare_abs(temp, remainder) <= 0) {
                left = mid;
            } else {
                right = mid - 1;
            }
        }
        
        quotient.digits.insert(quotient.digits.begin(), left);
        
        // Subtract left * divisor from remainder
        int2048 temp = divisor;
        long long carry = 0;
        for (size_t j = 0; j < temp.digits.size(); ++j) {
            long long prod = (long long)temp.digits[j] * left + carry;
            temp.digits[j] = prod % BASE;
            carry = prod / BASE;
        }
        if (carry > 0) {
            temp.digits.push_back(carry);
        }
        remainder = sub_abs(remainder, temp);
    }
    
    *this = quotient;
    remove_leading_zeros();
    
    // Handle floor division for negative results
    if (result_negative) {
        if (!(remainder.digits.size() == 1 && remainder.digits[0] == 0)) {
            // Has remainder, need to round down (more negative)
            *this += int2048(1);
        }
        this->negative = true;
    }
    
    return *this;
}

int2048 operator/(int2048 a, const int2048 &b) {
    a /= b;
    return a;
}

// Modulo
int2048 &int2048::operator%=(const int2048 &other) {
    int2048 quotient = *this / other;
    *this -= quotient * other;
    return *this;
}

int2048 operator%(int2048 a, const int2048 &b) {
    a %= b;
    return a;
}

// I/O operators
std::istream &operator>>(std::istream &is, int2048 &num) {
    std::string s;
    is >> s;
    num.read(s);
    return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &num) {
    if (num.negative) os << '-';
    os << num.digits.back();
    for (int i = num.digits.size() - 2; i >= 0; --i) {
        int d = num.digits[i];
        for (int j = BASE_DIGITS - 1; j > 0; --j) {
            if (d < (int)pow(10, j)) os << '0';
            else break;
        }
        os << d;
    }
    return os;
}

// Comparison operators
bool operator==(const int2048 &a, const int2048 &b) {
    return a.negative == b.negative && a.digits == b.digits;
}

bool operator!=(const int2048 &a, const int2048 &b) {
    return !(a == b);
}

bool operator<(const int2048 &a, const int2048 &b) {
    if (a.negative != b.negative) {
        return a.negative;
    }
    int cmp = int2048::compare_abs(a, b);
    return a.negative ? (cmp > 0) : (cmp < 0);
}

bool operator>(const int2048 &a, const int2048 &b) {
    return b < a;
}

bool operator<=(const int2048 &a, const int2048 &b) {
    return !(a > b);
}

bool operator>=(const int2048 &a, const int2048 &b) {
    return !(a < b);
}

} // namespace sjtu

#endif
